"""Command line interface for the VariantGrid API client.

Currently provides `vg_api annotate_vcf <file>`: a stateless upload/download helper. Because the server
dedups uploads on the file's SHA-256, the VCF file itself is the receipt - the first call uploads it, and
running the same command again downloads the annotated result (or reports that it isn't ready yet). No local
state (upload ids, pending.json, ...) is kept, so you can poll from any machine that has the VCF.
"""
import argparse
import hashlib
import json
import logging
import os
import sys

import requests

from variantgrid_api.api_client import VariantGridAPI, AnnotationError, UnsupportedFeatureError

DEFAULT_SERVER = "https://variantgrid.com"

EXIT_OK = 0        # annotated file downloaded
EXIT_ERROR = 1     # bad input / auth / server or annotation error
EXIT_PENDING = 3   # uploaded just now, or still annotating - come back later


def _sha256(path):
    """Content hash the server dedups on - identical to `sha256sum <file>`."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _enable_verbose_logging():
    """Send the vg_api loggers (the CLI's and the client's requests/responses) to stderr."""
    vg_api_logger = logging.getLogger("vg_api")
    vg_api_logger.setLevel(logging.DEBUG)
    # main() can run more than once in a process (tests) - reuse our handler, pointed at the current stderr
    for handler in vg_api_logger.handlers:
        if getattr(handler, "_vg_api_verbose", False):
            handler.setStream(sys.stderr)
            return
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    handler._vg_api_verbose = True
    vg_api_logger.addHandler(handler)


def _build_api(args):
    token = args.token or os.environ.get("VARIANTGRID_API_TOKEN")
    if not token:
        raise SystemExit("No API token - pass --token or set VARIANTGRID_API_TOKEN")
    server = args.server or os.environ.get("VARIANTGRID_API_SERVER") or DEFAULT_SERVER
    # Own logger so we can silence the expected 404 when a file hasn't been uploaded yet.
    logger = logging.getLogger("vg_api.cli")
    verbose = args.verbose
    if verbose:
        _enable_verbose_logging()
    api = VariantGridAPI(server=server, api_token=token, logger=logger,
                         log_request=verbose, log_response=verbose)
    return api, logger


def _log_status(args, logger, status):
    if args.verbose:
        logger.debug("Upload status: %s", json.dumps(status, indent=2, sort_keys=True))


def _pending_reason(status):
    """Which of the gates before download hasn't passed yet, from the upload_status dict, eg
    'pipeline Success, import 100%, 2 annotation runs remaining'."""
    parts = []
    if pipeline_status := status.get("pipeline_status"):
        parts.append(f"pipeline {pipeline_status}")
    import_status = status.get("import_status")
    progress = status.get("progress_percent")
    if import_status or progress is not None:
        import_part = "import"
        if import_status:
            import_part += f" {import_status}"
        if progress is not None:
            import_part += f" {progress}%"
        parts.append(import_part)
    remaining = status.get("remaining_annotation_runs")
    if remaining:
        parts.append(f"{remaining} annotation run{'s' if remaining != 1 else ''} remaining")
    elif status.get("annotation_complete") is False:
        parts.append("annotation not complete")
    if status.get("downloads_available") is False:
        parts.append("downloads not available")
    return ", ".join(parts)


def _is_not_found(exc):
    resp = getattr(exc, "response", None)
    return resp is not None and resp.status_code == 404


def _probe_status(api, logger, sha256, verbose=False):
    """Poll status by content hash. Returns the status dict, or None if the server has never
    seen this file (404) - a 404 here just means "not uploaded yet", so we silence its logging
    unless verbose."""
    prev_level = logger.level
    if not verbose:
        logger.setLevel(logging.CRITICAL)
    try:
        return api.poll_upload_status(sha256=sha256)
    except requests.HTTPError as e:
        if _is_not_found(e):
            return None
        raise
    finally:
        logger.setLevel(prev_level)


def _download(api, args, sha256):
    path = api.download_annotated(sha256=sha256, export_type=args.export_type, dest_path=args.dest)
    print(f"Annotated {args.export_type} written to {path}")
    return EXIT_OK


def annotate_vcf_cmd(args):
    if not os.path.isfile(args.vcf):
        print(f"No such file: {args.vcf}", file=sys.stderr)
        return EXIT_ERROR

    api, logger = _build_api(args)
    name = os.path.basename(args.vcf)
    sha256 = _sha256(args.vcf)

    status = _probe_status(api, logger, sha256, verbose=args.verbose)
    if status is None:
        # We've never uploaded this file - do it now.
        up = api.upload_file(args.vcf, path=None)
        print(f"Uploaded {name} (id={up['uploaded_file_id']}).")

    if args.wait:
        # Poll the server ourselves until annotation finishes (can take hours), then download.
        if status is None or not status.get("annotation_complete"):
            print("Waiting for annotation to finish - this can take a while...")
            api.wait_for_annotation(sha256=sha256, timeout=args.timeout, poll_interval=args.poll_interval)
        return _download(api, args, sha256)

    # Single-shot: report where it's at, download only if it's ready.
    if status is None:
        print(f"Annotating {name} - run the same command again later to download.")
        return EXIT_PENDING
    if err := status.get("error"):
        _log_status(args, logger, status)
        print(f"{name}: annotation error - {err}", file=sys.stderr)
        return EXIT_ERROR
    if status.get("annotation_complete"):
        return _download(api, args, sha256)

    _log_status(args, logger, status)
    reason = _pending_reason(status)
    reason = f" - {reason}" if reason else ""
    print(f"{name}: not ready yet{reason}. Run the same command again later.")
    return EXIT_PENDING


def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--server",
                        help="VariantGrid server URL (default: $VARIANTGRID_API_SERVER or "
                             f"{DEFAULT_SERVER})")
    common.add_argument("--token", help="API token (default: $VARIANTGRID_API_TOKEN)")
    common.add_argument("-v", "--verbose", action="store_true",
                        help="Log requests, responses and the full upload status to stderr")

    parser = argparse.ArgumentParser(prog="vg_api", description="VariantGrid API command line tool")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("annotate_vcf", parents=[common],
                       help="Upload a VCF for annotation; run again to download it once ready",
                       description="Upload a VCF for annotation. Because uploads are keyed on the file's "
                                   "SHA-256, running the same command again downloads the annotated result "
                                   "when ready, or reports that it isn't ready yet - no local state is kept.")
    p.add_argument("vcf", help="Path to the VCF (.vcf / .vcf.gz) to annotate")
    p.add_argument("--export-type", choices=("vcf", "csv"), default="vcf",
                   help="Download format (default: vcf)")
    p.add_argument("-o", "--dest", default=".",
                   help="Destination directory or file for the download (default: current directory)")
    p.add_argument("--wait", action="store_true",
                   help="Poll until annotation finishes and download it, instead of returning immediately "
                        "(can take hours)")
    p.add_argument("--poll-interval", type=float, default=10,
                   help="Seconds between status polls when using --wait (default: 10)")
    p.add_argument("--timeout", type=float, default=86400,
                   help="Give up after this many seconds when using --wait (default: 86400 = 24h)")
    p.set_defaults(func=annotate_vcf_cmd)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except AnnotationError as e:
        if e.status is not None:
            _log_status(args, logging.getLogger("vg_api.cli"), e.status)
        print(f"Annotation error: {e}", file=sys.stderr)
        return EXIT_ERROR
    except requests.HTTPError as e:
        print(f"HTTP error: {e}", file=sys.stderr)
        return EXIT_ERROR
    except UnsupportedFeatureError as e:
        print(f"This VariantGrid server can't annotate uploaded VCFs via the API: {e}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
