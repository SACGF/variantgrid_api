"""Command line interface for the VariantGrid API client.

Currently provides `vg_api annotate_vcf <file>`: a stateless upload/download helper. Because the server
dedups uploads on the file's SHA-256, the VCF file itself is the receipt - the first call uploads it, and
running the same command again downloads the annotated result (or reports that it isn't ready yet). No local
state (upload ids, pending.json, ...) is kept, so you can poll from any machine that has the VCF.
"""
import argparse
import hashlib
import logging
import os
import sys

import requests

from variantgrid_api.api_client import VariantGridAPI, AnnotationError

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


def _build_api(args):
    token = args.token or os.environ.get("VARIANTGRID_API_TOKEN")
    if not token:
        raise SystemExit("No API token - pass --token or set VARIANTGRID_API_TOKEN")
    server = args.server or os.environ.get("VARIANTGRID_API_SERVER") or DEFAULT_SERVER
    # Own logger so we can silence the expected 404 when a file hasn't been uploaded yet.
    logger = logging.getLogger("vg_api.cli")
    return VariantGridAPI(server=server, api_token=token, logger=logger), logger


def _is_not_found(exc):
    resp = getattr(exc, "response", None)
    return resp is not None and resp.status_code == 404


def annotate_vcf_cmd(args):
    if not os.path.isfile(args.vcf):
        print(f"No such file: {args.vcf}", file=sys.stderr)
        return EXIT_ERROR

    api, logger = _build_api(args)
    name = os.path.basename(args.vcf)

    if args.wait:
        # Blocking one-shot: upload, wait (can take hours), download.
        path = api.annotate_vcf(args.vcf, export_type=args.export_type, dest_path=args.dest,
                                poll_interval=args.poll_interval)
        print(f"Annotated {args.export_type} written to {path}")
        return EXIT_OK

    sha256 = _sha256(args.vcf)

    # Probe by content hash. A 404 means we've never uploaded this file - so upload it now.
    prev_level = logger.level
    logger.setLevel(logging.CRITICAL)  # the probe 404 is expected; don't scare the user
    try:
        status = api.poll_upload_status(sha256=sha256)
    except requests.HTTPError as e:
        if _is_not_found(e):
            up = api.upload_file(args.vcf, path=None)
            print(f"Uploaded {name} (id={up['uploaded_file_id']}). "
                  f"Annotating - run the same command again later to download.")
            return EXIT_PENDING
        raise
    finally:
        logger.setLevel(prev_level)

    if err := status.get("error"):
        print(f"{name}: annotation error - {err}", file=sys.stderr)
        return EXIT_ERROR
    if status.get("annotation_complete"):
        path = api.download_annotated(sha256=sha256, export_type=args.export_type, dest_path=args.dest)
        print(f"Annotated {args.export_type} written to {path}")
        return EXIT_OK

    progress = status.get("progress_percent")
    suffix = f" (progress {progress}%)" if progress is not None else ""
    print(f"{name}: not ready yet{suffix} - run the same command again later.")
    return EXIT_PENDING


def build_parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--server",
                        help="VariantGrid server URL (default: $VARIANTGRID_API_SERVER or "
                             f"{DEFAULT_SERVER})")
    common.add_argument("--token", help="API token (default: $VARIANTGRID_API_TOKEN)")

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
                   help="Block until annotation finishes, then download (can take hours)")
    p.add_argument("--poll-interval", type=float, default=10,
                   help="Seconds between status polls when using --wait (default: 10)")
    p.set_defaults(func=annotate_vcf_cmd)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except AnnotationError as e:
        print(f"Annotation error: {e}", file=sys.stderr)
        return EXIT_ERROR
    except requests.HTTPError as e:
        print(f"HTTP error: {e}", file=sys.stderr)
        return EXIT_ERROR


if __name__ == "__main__":
    sys.exit(main())
