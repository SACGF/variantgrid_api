import gzip

import pytest
import responses

from variantgrid_api.api_client import AnnotationError, VariantGridAPI


def _no_sleep(_seconds):
    """Injectable sleep that never actually waits."""
    return None


@responses.activate
def test_upload_file_omits_path_when_none(api, server, tmp_path):
    """Ad-hoc uploads must omit the SeqAuto-only 'path' query param."""
    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")
    url = f"{server}/upload/api/v1/file_upload"
    responses.add(responses.POST, url, json={"uploaded_file_id": 1, "sha256_hash": "abc"}, status=200)

    out = api.upload_file(str(src), path=None)
    assert out["uploaded_file_id"] == 1
    assert "path=" not in responses.calls[-1].request.url


@responses.activate
def test_upload_file_defaults_to_sending_path(api, server, tmp_path):
    """Backwards-compatible default (SeqAuto) still sends path=filename."""
    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")
    url = f"{server}/upload/api/v1/file_upload"
    responses.add(responses.POST, url, json={"uploaded_file_id": 1}, status=200)

    api.upload_file(str(src))
    assert "path=" in responses.calls[-1].request.url


@responses.activate
def test_annotate_vcf_uploads_without_path(api, server, tmp_path):
    """The convenience wrapper must upload ad-hoc (no path) so SeqAuto deployments don't reject it."""
    upload_url = f"{server}/upload/api/v1/file_upload"
    status_url = f"{server}/upload/api/v1/upload_status/789"
    download_url = f"{server}/upload/api/v1/download/789/vcf"
    responses.add(responses.POST, upload_url, json={"uploaded_file_id": 789}, status=200)
    responses.add(responses.GET, status_url, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, download_url, body=b"x", status=200,
                  headers={"Content-Disposition": 'attachment; filename="a.vcf.gz"'})

    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")
    api.annotate_vcf(str(src), dest_path=tmp_path, poll_interval=0, sleep=_no_sleep)

    upload_call = next(c for c in responses.calls if c.request.url.startswith(upload_url))
    assert "path=" not in upload_call.request.url


@responses.activate
def test_poll_upload_status_by_id(api, server):
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"uploaded_file_id": 123, "annotation_complete": False}, status=200)
    out = api.poll_upload_status(uploaded_file_id=123)
    assert out["uploaded_file_id"] == 123
    assert responses.calls[-1].request.headers["Authorization"] == "Token TKN"


@responses.activate
def test_poll_upload_status_by_sha256(api, server):
    url = f"{server}/upload/api/v1/upload_status/sha256/abc123"
    responses.add(responses.GET, url, json={"sha256_hash": "abc123", "annotation_complete": True}, status=200)
    out = api.poll_upload_status(sha256="abc123")
    assert out["annotation_complete"] is True
    assert responses.calls[-1].request.url.endswith("/upload_status/sha256/abc123")


def test_poll_upload_status_requires_a_key(api):
    with pytest.raises(ValueError):
        api.poll_upload_status()


@responses.activate
def test_wait_for_annotation_polls_until_complete(api, server):
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"annotation_complete": False, "error": None}, status=200)
    responses.add(responses.GET, url, json={"annotation_complete": False, "error": None}, status=200)
    responses.add(responses.GET, url, json={"annotation_complete": True, "error": None}, status=200)

    out = api.wait_for_annotation(uploaded_file_id=123, poll_interval=0, sleep=_no_sleep)
    assert out["annotation_complete"] is True
    assert len(responses.calls) == 3


@responses.activate
def test_wait_for_annotation_tolerates_transient_5xx(api, server):
    """A brief 500 (e.g. the upload->record race right after upload) should be retried, not fatal."""
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"detail": "server error"}, status=500)
    responses.add(responses.GET, url, json={"annotation_complete": False, "error": None}, status=200)
    responses.add(responses.GET, url, json={"annotation_complete": True, "error": None}, status=200)

    out = api.wait_for_annotation(uploaded_file_id=123, poll_interval=0, sleep=_no_sleep)
    assert out["annotation_complete"] is True


@responses.activate
def test_wait_for_annotation_gives_up_after_repeated_5xx(api, server):
    url = f"{server}/upload/api/v1/upload_status/123"
    for _ in range(6):
        responses.add(responses.GET, url, json={"detail": "server error"}, status=500)

    with pytest.raises(Exception):
        api.wait_for_annotation(uploaded_file_id=123, poll_interval=0, sleep=_no_sleep, max_transient_errors=3)


@responses.activate
def test_wait_for_annotation_does_not_retry_4xx(api, server):
    """A 4xx is a real error (bad key, not found) and should surface immediately."""
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"detail": "not found"}, status=404)

    with pytest.raises(Exception):
        api.wait_for_annotation(uploaded_file_id=123, poll_interval=0, sleep=_no_sleep)
    assert len(responses.calls) == 1  # no retry


@responses.activate
def test_wait_for_annotation_raises_on_error(api, server):
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"annotation_complete": False, "error": "Pipeline failed"}, status=200)

    with pytest.raises(AnnotationError) as exc_info:
        api.wait_for_annotation(uploaded_file_id=123, sleep=_no_sleep)
    assert "Pipeline failed" in str(exc_info.value)
    assert exc_info.value.status["error"] == "Pipeline failed"


@responses.activate
def test_wait_for_annotation_times_out(api, server):
    url = f"{server}/upload/api/v1/upload_status/123"
    responses.add(responses.GET, url, json={"annotation_complete": False, "error": None}, status=200)

    with pytest.raises(TimeoutError):
        api.wait_for_annotation(uploaded_file_id=123, timeout=0, sleep=_no_sleep)


@responses.activate
def test_download_annotated_streams_to_dest(api, server, tmp_path):
    url = f"{server}/upload/api/v1/download/123/vcf"
    payload = gzip.compress(b"##fileformat=VCFv4.2\n")
    responses.add(responses.GET, url, body=payload, status=200,
                  content_type="application/gzip",
                  headers={"Content-Disposition": 'attachment; filename="annotated.vcf.gz"'})

    dest = tmp_path / "out.vcf.gz"
    out = api.download_annotated(uploaded_file_id=123, export_type="vcf", dest_path=dest, sleep=_no_sleep)
    assert out == dest
    assert dest.read_bytes() == payload


@responses.activate
def test_download_annotated_uses_attachment_name_in_dir(api, server, tmp_path):
    url = f"{server}/upload/api/v1/download/123/csv"
    responses.add(responses.GET, url, body=b"zipbytes", status=200,
                  headers={"Content-Disposition": 'attachment; filename="cohort.csv.zip"'})

    out = api.download_annotated(uploaded_file_id=123, export_type="csv", dest_path=tmp_path, sleep=_no_sleep)
    assert out == tmp_path / "cohort.csv.zip"
    assert out.read_bytes() == b"zipbytes"


@responses.activate
def test_download_annotated_follows_202_generating(api, server, tmp_path):
    url = f"{server}/upload/api/v1/download/123/vcf"
    responses.add(responses.GET, url, json={"status": "generating", "progress": 0.5}, status=202)
    responses.add(responses.GET, url, body=b"data", status=200,
                  headers={"Content-Disposition": 'attachment; filename="a.vcf.gz"'})

    out = api.download_annotated(uploaded_file_id=123, dest_path=tmp_path, poll_interval=0, sleep=_no_sleep)
    assert out.read_bytes() == b"data"
    assert len(responses.calls) == 2


@responses.activate
def test_download_annotated_times_out_while_generating(api, server, tmp_path):
    url = f"{server}/upload/api/v1/download/123/vcf"
    responses.add(responses.GET, url, json={"status": "generating", "progress": 0.1}, status=202)

    with pytest.raises(TimeoutError):
        api.download_annotated(uploaded_file_id=123, dest_path=tmp_path, timeout=0, sleep=_no_sleep)


@responses.activate
def test_download_annotated_raises_on_server_error(api, server, tmp_path):
    url = f"{server}/upload/api/v1/download/123/vcf"
    responses.add(responses.GET, url, json={"error": "Export templates not configured"}, status=400)

    with pytest.raises(Exception):
        api.download_annotated(uploaded_file_id=123, dest_path=tmp_path, sleep=_no_sleep)


def test_download_annotated_rejects_bad_export_type(api):
    with pytest.raises(ValueError):
        api.download_annotated(uploaded_file_id=123, export_type="bam")


@responses.activate
def test_annotate_vcf_chains_upload_wait_download(api, server, tmp_path):
    upload_url = f"{server}/upload/api/v1/file_upload"
    status_url = f"{server}/upload/api/v1/upload_status/456"
    download_url = f"{server}/upload/api/v1/download/456/vcf"
    responses.add(responses.POST, upload_url, json={"uploaded_file_id": 456}, status=200)
    responses.add(responses.GET, status_url, json={"annotation_complete": True, "error": None}, status=200)
    responses.add(responses.GET, download_url, body=b"vcfdata", status=200,
                  headers={"Content-Disposition": 'attachment; filename="final.vcf.gz"'})

    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")

    out = api.annotate_vcf(str(src), dest_path=tmp_path, poll_interval=0, sleep=_no_sleep)
    assert out == tmp_path / "final.vcf.gz"
    assert out.read_bytes() == b"vcfdata"
