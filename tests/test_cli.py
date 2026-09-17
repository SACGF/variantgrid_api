import re

import pytest
import responses

from variantgrid_api import cli

SERVER = "https://vg.example.org"

STATUS_RE = re.compile(rf"{SERVER}/upload/api/v1/upload_status/sha256/[0-9a-f]+$")
UPLOAD_URL = f"{SERVER}/upload/api/v1/file_upload"
DOWNLOAD_RE = re.compile(rf"{SERVER}/upload/api/v1/download/sha256/[0-9a-f]+/(vcf|csv)$")
CAPABILITIES_URL = f"{SERVER}/seqauto/api/v1/capabilities"


def _add_capabilities(features=("upload_status",)):
    responses.add(responses.GET, CAPABILITIES_URL, status=200,
                  json={"version": "4.0.0", "git_hash": None, "features": list(features), "upload_file_types": ["vcf"]})


@pytest.fixture
def vcf(tmp_path):
    p = tmp_path / "input.vcf"
    p.write_text("##fileformat=VCFv4.2\n")
    return str(p)


def _argv(vcf, *extra):
    return ["annotate_vcf", vcf, "--server", SERVER, "--token", "T", *extra]


@responses.activate
def test_first_call_uploads(vcf):
    """Unknown file (status 404) -> upload it and report pending."""
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE, json={"detail": "not found"}, status=404)
    responses.add(responses.POST, UPLOAD_URL, json={"uploaded_file_id": 42}, status=200)

    rc = cli.main(_argv(vcf))

    assert rc == cli.EXIT_PENDING
    assert any(c.request.method == "POST" and c.request.url.startswith(UPLOAD_URL)
               for c in responses.calls)
    # Ad-hoc upload must omit the SeqAuto path hint
    upload_call = next(c for c in responses.calls if c.request.url.startswith(UPLOAD_URL))
    assert "path=" not in upload_call.request.url


@responses.activate
def test_second_call_downloads_when_ready(vcf, tmp_path):
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, DOWNLOAD_RE, body=b"vcfdata", status=200,
                  headers={"Content-Disposition": 'attachment; filename="out.vcf.gz"'})

    rc = cli.main(_argv(vcf, "-o", str(tmp_path)))

    assert rc == cli.EXIT_OK
    assert (tmp_path / "out.vcf.gz").read_bytes() == b"vcfdata"
    # It must NOT re-upload when the file is already known
    assert not any(c.request.method == "POST" for c in responses.calls)


@responses.activate
def test_not_ready_reports_pending(vcf):
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE,
                  json={"annotation_complete": False, "error": None, "progress_percent": 40}, status=200)

    rc = cli.main(_argv(vcf))

    assert rc == cli.EXIT_PENDING
    assert not any(re.match(DOWNLOAD_RE, c.request.url) for c in responses.calls)


@responses.activate
def test_server_error_status_is_error(vcf):
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE,
                  json={"annotation_complete": False, "error": "Pipeline failed"}, status=200)

    rc = cli.main(_argv(vcf))
    assert rc == cli.EXIT_ERROR


@responses.activate
def test_csv_export_type(vcf, tmp_path):
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, DOWNLOAD_RE, body=b"zip", status=200,
                  headers={"Content-Disposition": 'attachment; filename="out.csv.zip"'})

    rc = cli.main(_argv(vcf, "--export-type", "csv", "-o", str(tmp_path)))

    assert rc == cli.EXIT_OK
    assert (tmp_path / "out.csv.zip").read_bytes() == b"zip"
    download_call = next(c for c in responses.calls if re.match(DOWNLOAD_RE, c.request.url))
    assert download_call.request.url.endswith("/csv")


@responses.activate
def test_wait_uploads_polls_then_downloads(vcf, tmp_path):
    """--wait on a new file: upload, poll until complete, then download - no external loop."""
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE, json={"detail": "not found"}, status=404)  # probe
    responses.add(responses.POST, UPLOAD_URL, json={"uploaded_file_id": 7}, status=200)
    responses.add(responses.GET, STATUS_RE, json={"annotation_complete": False, "error": None}, status=200)
    responses.add(responses.GET, STATUS_RE, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, DOWNLOAD_RE, body=b"data", status=200,
                  headers={"Content-Disposition": 'attachment; filename="out.vcf.gz"'})

    rc = cli.main(_argv(vcf, "--wait", "-o", str(tmp_path), "--poll-interval", "0"))

    assert rc == cli.EXIT_OK
    assert (tmp_path / "out.vcf.gz").read_bytes() == b"data"


@responses.activate
def test_wait_downloads_without_reupload_when_already_known(vcf, tmp_path):
    """--wait on an already-uploaded, already-complete file must not re-upload."""
    _add_capabilities()
    responses.add(responses.GET, STATUS_RE, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, DOWNLOAD_RE, body=b"data", status=200,
                  headers={"Content-Disposition": 'attachment; filename="out.vcf.gz"'})

    rc = cli.main(_argv(vcf, "--wait", "-o", str(tmp_path), "--poll-interval", "0"))

    assert rc == cli.EXIT_OK
    assert not any(c.request.method == "POST" for c in responses.calls)


def test_missing_file_is_error(tmp_path):
    rc = cli.main(_argv(str(tmp_path / "nope.vcf")))
    assert rc == cli.EXIT_ERROR


def test_missing_token_raises(vcf, monkeypatch):
    monkeypatch.delenv("VARIANTGRID_API_TOKEN", raising=False)
    with pytest.raises(SystemExit):
        cli.main(["annotate_vcf", vcf, "--server", SERVER])


@responses.activate
def test_server_without_upload_status_is_error(vcf, capsys):
    """ A server that predates the capabilities endpoint (404) can't be used - nothing is uploaded """
    responses.add(responses.GET, CAPABILITIES_URL, status=404)

    rc = cli.main(_argv(vcf))

    assert rc == cli.EXIT_ERROR
    assert "can't annotate" in capsys.readouterr().err
    assert len(responses.calls) == 1
