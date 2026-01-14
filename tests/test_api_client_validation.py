import pytest, responses
from variantgrid_api.api_client import VariantGridAPI, EmptyInputPolicy

import responses

def run_validation_call(api_call, url=None, method=responses.POST, response_json={"ok": True}):
    if url:
        responses.add(method, url, json=response_json, status=200)
    return api_call()

@responses.activate
def test_error_policy_empty_list_raises(server, api_token):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.ERROR)

    with pytest.raises(ValueError):
        run_validation_call(
            lambda: api.create_multiple_qc_gene_lists([])
        )

    assert len(responses.calls) == 0

@responses.activate
def test_warn_policy_empty_list_logs_and_posts(server, api_token, caplog):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.WARN)
    url = f"{server}/seqauto/api/v1/qc_gene_list/bulk_create"

    out = run_validation_call(
        lambda: api.create_multiple_qc_gene_lists([]),
        url=url
    )

    assert out == {"ok": True}
    assert len(responses.calls) == 1
    assert any("empty list" in r.message for r in caplog.records)

@responses.activate
def test_ignore_policy_empty_list_posts_silently(server, api_token, caplog):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.IGNORE)
    url = f"{server}/seqauto/api/v1/qc_gene_list/bulk_create"

    out = run_validation_call(
        lambda: api.create_multiple_qc_gene_lists([]),
        url=url
    )

    assert out == {"ok": True}
    assert len(responses.calls) == 1
    assert not caplog.records


import pytest, responses
from variantgrid_api.api_client import VariantGridAPI, EmptyInputPolicy

@responses.activate
def test_validate_string_empty_raises(server, api_token):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.ERROR)
    with pytest.raises(ValueError):
        api.create_experiment("")
    assert len(responses.calls) == 0


@responses.activate
def test_validate_string_none_raises(server, api_token):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.ERROR)
    with pytest.raises(ValueError):
        api.create_experiment(None)
    assert len(responses.calls) == 0


@responses.activate
def test_validate_object_none_raises(server, api_token):
    api = VariantGridAPI(server, api_token, empty_input_policy=EmptyInputPolicy.ERROR)
    with pytest.raises(ValueError):
        api.create_enrichment_kit(None)
    assert len(responses.calls) == 0
