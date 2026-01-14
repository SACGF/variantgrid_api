import datetime
import json
import logging

import pytest
import responses

from variantgrid_api.api_client import VariantGridAPI, DateTimeEncoder


def _last_json():
    r = responses.calls[-1].request
    return json.loads(r.body)

@responses.activate
def test_create_experiment(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/experiment/"
    responses.add(responses.POST, url, json={"id": 1}, status=200)
    out = api.create_experiment(vg_objects["experiment"])
    assert out == {"id": 1}
    assert responses.calls[-1].request.headers["Authorization"] == "Token TKN"

@responses.activate
def test_create_sample_sheet_rewrites_sequencing_run_to_name(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sample_sheet/"
    responses.add(responses.POST, url, json={"id": 2}, status=200)
    api.create_sample_sheet(vg_objects["sample_sheet"])
    body = _last_json()
    assert isinstance(body["sequencing_run"], str)

@responses.activate
def test_create_sequencing_data_builds_unaligned_reads(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sequencing_files/bulk_create"
    responses.add(responses.POST, url, json={"created": 2}, status=200)
    api.create_sequencing_data(vg_objects["sample_sheet_lookup"], vg_objects["sequencing_files"])
    body = _last_json()
    assert "records" in body and len(body["records"]) == 2
    r0 = body["records"][0]
    assert "unaligned_reads" in r0
    assert "fastq_r1" in r0["unaligned_reads"] and "path" in r0["unaligned_reads"]["fastq_r1"]

def assert_post(api_call, url):
    responses.add(responses.POST, url, json={"ok": True}, status=200)
    out = api_call()
    assert out == {"ok": True}
    body = json.loads(responses.calls[-1].request.body)
    assert body
    return body

@responses.activate
def test_create_enrichment_kit_posts_json(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/enrichment_kit/"
    body = assert_post(
        lambda: api.create_enrichment_kit(vg_objects["enrichment_kit"]),
        url
    )
    assert body["name"] == "idt_haem"


@responses.activate
def test_create_sequencer_model_posts_json(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sequencer_model/"
    body = assert_post(
        lambda: api.create_sequencer_model(vg_objects["sequencer_model"]),
        url
    )
    assert body["model"] == "HiSeq 2500"


@responses.activate
def test_create_sequencer_posts_json(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sequencer/"
    body = assert_post(
        lambda: api.create_sequencer(vg_objects["sequencer"]),
        url
    )
    assert body["name"] == "SN1101"


@responses.activate
def test_create_sequencing_run_posts_json(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sequencing_run/"
    body = assert_post(
        lambda: api.create_sequencing_run(vg_objects["sequencing_run"]),
        url
    )
    assert body["name"] == vg_objects["SEQUENCING_RUN_NAME"]


@responses.activate
def test_create_sample_sheet_combined_vcf_file_posts_json(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sample_sheet_combined_vcf_file/"
    body = assert_post(
        lambda: api.create_sample_sheet_combined_vcf_file(
            vg_objects["sample_sheet_combined_vcf_file"]
        ),
        url
    )
    assert body["path"].endswith(".vcf.gz")


@responses.activate
def test_upload_file_calls_expected_url_and_query_param(api, server, vg_objects):
    filename = vg_objects["qc_gene_coverage_list"][0].path
    url = f"{server}/upload/api/v1/file_upload"
    responses.add(responses.POST, url, json={"uploaded": True}, status=200)
    out = api.upload_file(filename)
    assert out == {"uploaded": True}
    assert responses.calls[-1].request.url.startswith(url)
    assert "path=" in responses.calls[-1].request.url


@responses.activate
def test_sequencing_run_has_vcf(api, server, vg_objects):
    sequencing_run = vg_objects["sequencing_run"]
    name = sequencing_run.name
    url = f"{server}/seqauto/api/v1/sequencing_run/{name}/"
    responses.add(responses.GET, url, json={"vcf_set":[{"path":"/x/a.vcf.gz"}]}, status=200)

    assert api.sequencing_run_has_vcf(sequencing_run) is True
    assert api.sequencing_run_has_vcf(sequencing_run, "a.vcf.gz") is True
    assert api.sequencing_run_has_vcf(sequencing_run, "no_vcf_here") is False


@responses.activate
def test_create_qc_gene_list_posts(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/qc_gene_list/"
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    out = api.create_qc_gene_list(vg_objects["qc_gene_lists"][0])
    assert out == {"ok": True}
    body = _last_json()
    assert "gene_list" in body and len(body["gene_list"]) > 0

@responses.activate
def test_create_qc_exec_stats_posts(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/qc_exec_summary/"
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    out = api.create_qc_exec_stats(vg_objects["qc_exec_stats"][0])
    assert out == {"ok": True}
    body = _last_json()
    assert "reads" in body

import responses

@responses.activate
def test_log_request_emits_log(server, api_token, vg_objects, caplog):
    logger = logging.getLogger("variantgrid_api.test")
    caplog.set_level(logging.INFO, logger=logger.name)

    api = VariantGridAPI(
        server, api_token,
        logger=logger,
        log_request=True
    )
    url = f"{server}/seqauto/api/v1/experiment/"
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    api.create_experiment(vg_objects["experiment"])

    assert any("POST to" in r.message for r in caplog.records)



@responses.activate
def test_log_response_emits_log(server, api_token, vg_objects, caplog):
    logger = logging.getLogger("variantgrid_api.test")
    caplog.set_level(logging.INFO, logger=logger.name)

    api = VariantGridAPI(
        server, api_token,
        logger=logger,
        log_response=True
    )
    url = f"{server}/seqauto/api/v1/experiment/"
    responses.add(responses.POST, url, json={"ok": True}, status=200)

    api.create_experiment(vg_objects["experiment"])

    assert any("Response from" in r.message for r in caplog.records)


@responses.activate
def test_bad_response_raises_and_logs(server, api_token, vg_objects, caplog):
    logger = logging.getLogger("variantgrid_api.test")
    caplog.set_level(logging.INFO, logger=logger.name)

    api = VariantGridAPI(server, api_token, logger=logger)
    url = f"{server}/seqauto/api/v1/experiment/"
    responses.add(
        responses.POST,
        url,
        json={"detail": "bad request"},
        status=400
    )

    with pytest.raises(Exception):
        api.create_experiment(vg_objects["experiment"])

    assert any("Response:" in r.message for r in caplog.records)


@responses.activate
def test_bad_response_non_json_body(server, api_token, vg_objects):
    api = VariantGridAPI(server, api_token)
    url = f"{server}/seqauto/api/v1/experiment/"
    responses.add(
        responses.POST,
        url,
        body="not json",
        status=500,
        content_type="text/plain"
    )

    with pytest.raises(Exception):
        api.create_experiment(vg_objects["experiment"])


def test_datetime_encoder_handles_date_and_rejects_unknown():
    s = json.dumps(
        {"d": datetime.date(2026, 1, 14)},
        cls=DateTimeEncoder
    )
    assert s == '{"d": "2026-01-14"}'

    with pytest.raises(TypeError):
        json.dumps({"x": object()}, cls=DateTimeEncoder)
