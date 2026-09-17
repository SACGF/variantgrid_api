"""Patient / specimen / extraction, specimen measures, extraction links and upload metadata
(SACGF/variantgrid_api#20)"""
import json
from urllib.parse import urlparse, parse_qs

import pytest
import responses

from variantgrid_api.data_models import ExternalReference, ExternalPK, Extraction, SpecimenMeasure, \
    SpecimenMeasureType


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
# Specimen measures                                                    #
# ------------------------------------------------------------------ #

@responses.activate
def test_create_specimen_measure(api, server, vg_objects):
    responses.add(responses.POST, f"{server}/patients/api/v1/specimen_measure/", json={"id": 4}, status=201)
    api.create_specimen_measure("2600000001", vg_objects["specimen_measures"][0])
    assert _last_json() == {
        "specimen": "2600000001",
        "measure_type": "T",
        "value": 7.1,
        "unit": "mut/Mb",
        "method": "DRAGEN TSO500 2.1.1",
        "source_payload": {"Total TMB": "7.1", "Coding Region Size in Megabases": "1.27"},
        "extraction": "2600000001C",
    }


@responses.activate
def test_create_specimen_measures_bulk(api, server, vg_objects):
    url = f"{server}/patients/api/v1/specimen_measure/bulk_create"
    responses.add(responses.POST, url, json={"specimen": "2600000001", "measures": ["x", "y"]}, status=201)
    specimen_reference = ExternalReference.from_specimen(vg_objects["specimen"])
    api.create_specimen_measures(specimen_reference, vg_objects["specimen_measures"])
    body = _last_json()
    assert body["specimen"] == "2600000001"
    assert [m["measure_type"] for m in body["measures"]] == ["T", "M"]
    msi = body["measures"][1]
    assert msi == {"measure_type": "M", "value": 2.48, "unit": "%", "call": "Stable",
                   "threshold": ">= 20.00%", "threshold_source": "Illumina"}


def test_create_specimen_measures_empty_list_raises(api):
    with pytest.raises(ValueError):
        api.create_specimen_measures("2600000001", [])


def test_create_specimen_measure_empty_reference_raises(api):
    with pytest.raises(ValueError):
        api.create_specimen_measure("", SpecimenMeasure(measure_type=SpecimenMeasureType.TMB, value=1.0))


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
