import json, responses

def assert_bulk_post(api_call, url, n):
    responses.add(responses.POST, url, json={"ok": True}, status=200)
    out = api_call()
    assert out == {"ok": True}
    body = json.loads(responses.calls[-1].request.body)
    assert "records" in body
    assert isinstance(body["records"], list)
    assert len(body["records"]) == n
    assert all(isinstance(r, dict) and r for r in body["records"])
    return body["records"]


@responses.activate
def test_create_multiple_qc_gene_lists(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/qc_gene_list/bulk_create"
    records = assert_bulk_post(
        lambda: api.create_multiple_qc_gene_lists(vg_objects["qc_gene_lists"]),
        url,
        n=len(vg_objects["qc_gene_lists"])
    )
    assert "gene_list" in records[0]
    assert len(records[0]["gene_list"]) > 0


@responses.activate
def test_create_multiple_qc_exec_stats(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/qc_exec_summary/bulk_create"
    records = assert_bulk_post(
        lambda: api.create_multiple_qc_exec_stats(vg_objects["qc_exec_stats"]),
        url,
        n=len(vg_objects["qc_exec_stats"])
    )
    assert "reads" in records[0]


@responses.activate
def test_create_multiple_qc_gene_coverage(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/qc_gene_coverage/bulk_create"
    records = assert_bulk_post(
        lambda: api.create_multiple_qc_gene_coverage(vg_objects["qc_gene_coverage_list"]),
        url,
        n=len(vg_objects["qc_gene_coverage_list"])
    )
    assert "path" in records[0]
