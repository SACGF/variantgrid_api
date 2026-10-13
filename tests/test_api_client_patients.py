"""Patient / specimen / extraction, specimen measures, extraction links and upload metadata
(SACGF/variantgrid_api#20)"""
import json
import logging
import re
from urllib.parse import urlparse, parse_qs

import pytest
import responses

from variantgrid_api.api_client import VariantGridAPI, UnsupportedFeaturePolicy, UnsupportedFeatureError
from variantgrid_api.data_models import ExternalReference, ExternalPK, Extraction, ServerCapabilities, \
    ServerFeature


@pytest.fixture
def api(vg4_api):
    return vg4_api


def _last_json():
    return json.loads(responses.calls[-1].request.body)


def _last_query_params():
    return {k: v[0] for k, v in parse_qs(urlparse(responses.calls[-1].request.url).query).items()}


# ------------------------------------------------------------------ #
# ExternalReference                                                    #
# ------------------------------------------------------------------ #

def test_external_reference_with_only_reference_id_is_bare_string():
    assert ExternalReference(reference_id="2600000001C").to_json_value() == "2600000001C"


def test_external_reference_with_code_is_object():
    ref = ExternalReference(code="H12345", external_type="HelixID")
    assert ref.to_json_value() == {"code": "H12345", "external_type": "HelixID"}


def test_external_reference_code_needs_external_type():
    with pytest.raises(ValueError):
        ExternalReference(code="H12345")


def test_external_reference_needs_reference_id_or_code():
    with pytest.raises(ValueError):
        ExternalReference(external_type="HelixID")


def test_external_reference_from_patient_uses_patient_code_and_external_pk(vg_objects):
    ref = ExternalReference.from_patient(vg_objects["patient"])
    assert ref.to_json_value() == {"reference_id": "C0000001", "code": "H12345",
                                   "external_type": "HelixID", "external_manager": "HELIX"}


def test_external_reference_from_extraction_without_external_pk(vg_objects):
    ref = ExternalReference.from_extraction(vg_objects["dna_extraction"])
    assert ref.to_json_value() == "2600000001C"


# ------------------------------------------------------------------ #
# Patient / specimen / extraction                                      #
# ------------------------------------------------------------------ #

@responses.activate
def test_create_patient(api, server, vg_objects):
    responses.add(responses.POST, f"{server}/patients/api/v1/patient/", json={"id": 1}, status=201)
    assert api.create_patient(vg_objects["patient"]) == {"id": 1}
    assert _last_json() == {
        "patient_code": "C0000001",
        "date_of_birth": "1970-01-01",
        "sex": "F",
        "affected": True,
        "external_pk": {"code": "H12345", "external_type": "HelixID", "external_manager": "HELIX"},
    }


@responses.activate
def test_create_specimen_names_patient_by_reference(api, server, vg_objects):
    responses.add(responses.POST, f"{server}/patients/api/v1/specimen/", json={"id": 2}, status=201)
    api.create_specimen(vg_objects["specimen"])
    body = _last_json()
    assert body["patient"] == {"reference_id": "C0000001", "code": "H12345",
                               "external_type": "HelixID", "external_manager": "HELIX"}
    assert body["reference_id"] == "2600000001"
    assert body["tissue_status"] == "A"
    assert body["collection_date"] == "2026-01-01T09:00:00+10:30"
    assert "description" not in body and "external_pk" not in body


@responses.activate
def test_create_extraction_bare_string_specimen(api, server, vg_objects):
    responses.add(responses.POST, f"{server}/patients/api/v1/extraction/", json={"id": 3}, status=201)
    api.create_extraction(vg_objects["dna_extraction"])
    assert _last_json() == {"specimen": "2600000001", "reference_id": "2600000001C", "nucleic_acid_source": "D"}


@responses.activate
def test_create_extraction_specimen_by_external_pk(api, server):
    responses.add(responses.POST, f"{server}/patients/api/v1/extraction/", json={"id": 3}, status=201)
    extraction = Extraction(specimen=ExternalReference(code="S1", external_type="LabNo"),
                            external_pk=ExternalPK(code="E1", external_type="LabNo", external_manager="HELIX"))
    api.create_extraction(extraction)
    assert _last_json() == {"specimen": {"code": "S1", "external_type": "LabNo"},
                            "external_pk": {"code": "E1", "external_type": "LabNo", "external_manager": "HELIX"}}


@responses.activate
def test_create_specimen_unknown_patient_raises(api, server, vg_objects):
    responses.add(responses.POST, f"{server}/patients/api/v1/specimen/",
                  json={"patient": ["No Patient found"]}, status=400)
    with pytest.raises(Exception):
        api.create_specimen(vg_objects["specimen"])


# ------------------------------------------------------------------ #
# Linking a sequencing sample to its extraction                        #
# ------------------------------------------------------------------ #

@responses.activate
def test_link_sequencing_sample_extraction(api, server, vg_objects):
    url = f"{server}/seqauto/api/v1/sequencing_sample/link_extraction"
    responses.add(responses.POST, url, json={"match_status": "Matched", "match_error": None}, status=200)
    out = api.link_sequencing_sample_extraction(vg_objects["sequencing_sample_lookup_1"], "2600000001C")
    assert out["match_status"] == "Matched"
    assert _last_json() == {
        "sequencing_sample": {
            "sample_sheet": vg_objects["sample_sheet_lookup"].to_dict(),
            "sample_name": "fake_sample_1",
        },
        "extraction": "2600000001C",
    }


@responses.activate
def test_link_sequencing_sample_extraction_pending_is_not_an_error(api, server, vg_objects):
    """ An extraction the server doesn't have yet is a 202 - the claim is parked, not rejected """
    url = f"{server}/seqauto/api/v1/sequencing_sample/link_extraction"
    pending = {"sequencing_sample": "fake_sample_1", "match_status": "Pending",
               "match_error": "No Extraction found for reference_id=2600000001C", "extraction": None}
    responses.add(responses.POST, url, json=pending, status=202)
    out = api.link_sequencing_sample_extraction(vg_objects["sequencing_sample_lookup_1"], "2600000001C")
    assert out == pending


# ------------------------------------------------------------------ #
# Linking a sample sheet's sequencing samples in one call (#32)         #
# ------------------------------------------------------------------ #

def _bulk_url(server):
    return f"{server}/seqauto/api/v1/sequencing_sample/link_extractions"


def _single_url(server):
    return f"{server}/seqauto/api/v1/sequencing_sample/link_extraction"


def _link_result(sample_name, extraction, match_status="Matched"):
    return {"sequencing_sample": sample_name, "match_status": match_status, "match_error": None,
            "extraction": extraction if match_status == "Matched" else None}


EXTRACTIONS = {"fake_sample_1": "2600000001C", "fake_sample_2": "2600000002C"}


@responses.activate
def test_link_sequencing_sample_extractions(api, server, vg_objects):
    result = [_link_result("fake_sample_1", "2600000001C"), _link_result("fake_sample_2", "2600000002C")]
    responses.add(responses.POST, _bulk_url(server), json=result, status=200)

    out = api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], EXTRACTIONS)
    assert out == result
    assert len(responses.calls) == 1
    assert _last_json() == {
        "sample_sheet": vg_objects["sample_sheet_lookup"].to_dict(),
        "links": [
            {"sample_name": "fake_sample_1", "extraction": "2600000001C"},
            {"sample_name": "fake_sample_2", "extraction": "2600000002C"},
        ],
    }


@responses.activate
def test_link_sequencing_sample_extractions_external_reference(api, server, vg_objects):
    """ Each link takes the same reference forms as the single call """
    responses.add(responses.POST, _bulk_url(server), json=[], status=200)
    lims = ExternalReference(code="LIMS-1", external_type="LabNo")
    api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], {"fake_sample_1": lims})
    assert _last_json()["links"] == [{"sample_name": "fake_sample_1", "extraction": lims.to_json_value()}]


@responses.activate
def test_link_sequencing_sample_extractions_pending_is_not_an_error(api, server, vg_objects):
    """ match_status is per sample - a 202 with one 'Pending' still links the others """
    result = [_link_result("fake_sample_1", "2600000001C"),
              _link_result("fake_sample_2", None, match_status="Pending")]
    responses.add(responses.POST, _bulk_url(server), json=result, status=202)
    assert api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], EXTRACTIONS) == result


@pytest.fixture
def link_extraction_only_api(server, api_token):
    """ A server from before the bulk endpoint, with just the single-sample one """
    api = VariantGridAPI(server, api_token)
    api._capabilities = ServerCapabilities(version="old", features=frozenset([ServerFeature.LINK_EXTRACTION]))
    return api


@responses.activate
def test_link_sequencing_sample_extractions_falls_back_to_one_call_per_sample(link_extraction_only_api, server,
                                                                              vg_objects):
    """ So pipelines can use the bulk call unconditionally - the results are collected in order """
    api = link_extraction_only_api
    responses.add(responses.POST, _single_url(server), json=_link_result("fake_sample_1", "2600000001C"), status=200)
    responses.add(responses.POST, _single_url(server), json=_link_result("fake_sample_2", None, "Pending"), status=202)

    out = api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], EXTRACTIONS)
    assert out == [_link_result("fake_sample_1", "2600000001C"), _link_result("fake_sample_2", None, "Pending")]
    assert [c.request.url for c in responses.calls] == [_single_url(server)] * 2
    sample_sheet = vg_objects["sample_sheet_lookup"].to_dict()
    assert [json.loads(c.request.body) for c in responses.calls] == [
        {"sequencing_sample": {"sample_sheet": sample_sheet, "sample_name": "fake_sample_1"},
         "extraction": "2600000001C"},
        {"sequencing_sample": {"sample_sheet": sample_sheet, "sample_name": "fake_sample_2"},
         "extraction": "2600000002C"},
    ]


@responses.activate
def test_link_sequencing_sample_extractions_policy_applies_only_without_link_extraction(server, api_token,
                                                                                        vg_objects, caplog):
    """ unsupported_feature_policy applies only when the server has neither feature """
    legacy = ServerCapabilities.LEGACY
    api = VariantGridAPI(server, api_token)
    api._capabilities = legacy
    with pytest.raises(UnsupportedFeatureError):
        api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], EXTRACTIONS)

    skip_api = VariantGridAPI(server, api_token, unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)
    skip_api._capabilities = legacy
    with caplog.at_level(logging.WARNING):
        assert skip_api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], EXTRACTIONS) is None
    assert "Skipping" in caplog.text
    assert not responses.calls


@pytest.mark.parametrize("extractions, message", [
    ({}, "extractions: empty"),
    ({"fake_sample_1": "2600000001C", "fake_sample_2": None}, "extractions['fake_sample_2']: is None"),
    ({"fake_sample_1": ""}, "extractions['fake_sample_1']: empty string"),
    ({"": "2600000001C"}, "extractions: sample name: empty string"),
])
@responses.activate
def test_link_sequencing_sample_extractions_validation_names_the_sample(api, vg_objects, extractions, message):
    with pytest.raises(ValueError, match=re.escape(message)):
        api.link_sequencing_sample_extractions(vg_objects["sample_sheet_lookup"], extractions)
    assert not responses.calls


@responses.activate
def test_link_sequencing_sample_extractions_requires_sample_sheet_lookup(api, vg_objects):
    with pytest.raises(ValueError, match="sample_sheet_lookup: is None"):
        api.link_sequencing_sample_extractions(None, EXTRACTIONS)


# ------------------------------------------------------------------ #
# Upload metadata                                                      #
# ------------------------------------------------------------------ #

@pytest.fixture
def vcf(tmp_path):
    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")
    return str(src)


@responses.activate
def test_upload_file_sends_metadata_as_query_params(api, server, vcf):
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)
    api.upload_file(vcf, path=None, metadata={"genome_build": "GRCh37",
                                              "source": "DRAGEN TSO500 SmallVariant",
                                              "extraction": "2600000001C"})
    assert _last_query_params() == {"genome_build": "GRCh37", "source": "DRAGEN TSO500 SmallVariant",
                                    "extraction": "2600000001C"}


@responses.activate
def test_upload_file_metadata_references_sent_as_json(api, server, vcf):
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)
    tumour = ExternalReference(code="E1", external_type="LabNo")
    api.upload_file(vcf, path=None, metadata={"sample_extractions": {"TUMOUR": tumour, "NORMAL": "2600000002C"}})
    params = _last_query_params()
    assert json.loads(params["sample_extractions"]) == {
        "TUMOUR": {"code": "E1", "external_type": "LabNo"},
        "NORMAL": "2600000002C",
    }


@responses.activate
def test_upload_file_metadata_external_reference_extraction(api, server, vcf):
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)
    api.upload_file(vcf, path=None, metadata={"extraction": ExternalReference(reference_id="2600000001C")})
    assert _last_query_params() == {"extraction": "2600000001C"}


@responses.activate
def test_upload_file_metadata_keeps_path(api, server, vcf):
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)
    api.upload_file(vcf, metadata={"genome_build": "GRCh37"})
    assert _last_query_params() == {"path": vcf, "genome_build": "GRCh37"}


@responses.activate
def test_upload_file_without_metadata_unchanged(api, server, vcf):
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)
    api.upload_file(vcf)
    assert _last_query_params() == {"path": vcf}


@pytest.mark.parametrize("reserved", ["path", "force"])
def test_upload_file_metadata_reserved_key_raises(api, vcf, reserved):
    with pytest.raises(ValueError):
        api.upload_file(vcf, metadata={reserved: "x"})
