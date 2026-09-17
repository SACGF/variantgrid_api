"""Server capabilities probe and feature-gated calls (SACGF/variantgrid_api#21, SACGF/variantgrid_sapath#443)"""
import logging

import pytest
import requests
import responses

from variantgrid_api.api_client import VariantGridAPI, UnsupportedFeaturePolicy, UnsupportedFeatureError
from variantgrid_api.data_models import ServerCapabilities

CVO = "dragen_tso500_combined_variant_output"


@pytest.fixture
def capabilities_url(server):
    return f"{server}/api/v1/capabilities"


@pytest.fixture
def skip_api(server, api_token):
    return VariantGridAPI(server, api_token, unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)


@pytest.fixture
def vcf(tmp_path):
    src = tmp_path / "input.vcf"
    src.write_text("##fileformat=VCFv4.2\n")
    return str(src)


def _capabilities_calls(url):
    return [c for c in responses.calls if c.request.url == url]


# ------------------------------------------------------------------ #
# Probe                                                                #
# ------------------------------------------------------------------ #

@responses.activate
def test_capabilities_parsed_and_fetched_once(api, capabilities_url, capabilities_json):
    responses.add(responses.GET, capabilities_url, json=capabilities_json, status=200)

    assert api.supports("patients")
    assert api.supports("upload_metadata")
    assert not api.supports("no_such_feature")
    assert api.accepts_upload(CVO)
    assert not api.accepts_upload("gene_list")

    capabilities = api.capabilities
    assert capabilities.version == "4.0.0"
    assert capabilities.git_hash == "2130cffe0"
    assert isinstance(capabilities.features, frozenset)
    assert len(_capabilities_calls(capabilities_url)) == 1


@responses.activate
def test_capabilities_404_is_legacy(api, capabilities_url):
    responses.add(responses.GET, capabilities_url, json={"detail": "Not found."}, status=404)

    assert api.capabilities is ServerCapabilities.LEGACY
    assert api.capabilities.version == "legacy"
    assert not api.supports("patients")
    assert not api.accepts_upload("vcf")
    assert len(_capabilities_calls(capabilities_url)) == 1


@responses.activate
def test_capabilities_login_redirect_is_legacy(api, server, capabilities_url):
    """ An older server whose login middleware doesn't exempt /api/ redirects to its login page """
    login_url = f"{server}/accounts/login/?next=/api/v1/capabilities"
    responses.add(responses.GET, capabilities_url, status=302, headers={"Location": login_url})
    responses.add(responses.GET, login_url, body="<html>login</html>", status=200)

    assert api.capabilities is ServerCapabilities.LEGACY
    assert len(responses.calls) == 1  # redirect not followed


@responses.activate
def test_capabilities_other_error_raises(api, capabilities_url):
    responses.add(responses.GET, capabilities_url, json={"detail": "Invalid token."}, status=401)

    with pytest.raises(requests.HTTPError):
        api.capabilities


def test_constructing_client_makes_no_request(server, api_token):
    with responses.RequestsMock() as rsps:  # raises if any request is made
        VariantGridAPI(server, api_token)
        assert len(rsps.calls) == 0


@responses.activate
def test_ungated_method_makes_no_capabilities_request(api, server, capabilities_url, vg_objects, vcf):
    responses.add(responses.POST, f"{server}/seqauto/api/v1/sequencing_run/", json={}, status=200)
    responses.add(responses.POST, f"{server}/upload/api/v1/file_upload", json={"uploaded_file_id": 1}, status=200)

    api.create_sequencing_run(vg_objects["sequencing_run"])
    api.upload_file(vcf)
    assert not _capabilities_calls(capabilities_url)


# ------------------------------------------------------------------ #
# Gated methods on a legacy server                                     #
# ------------------------------------------------------------------ #

def _gated_calls(vg_objects, vcf):
    measure = vg_objects["specimen_measures"][0]
    return {
        "create_patient": lambda api: api.create_patient(vg_objects["patient"]),
        "create_specimen": lambda api: api.create_specimen(vg_objects["specimen"]),
        "create_extraction": lambda api: api.create_extraction(vg_objects["dna_extraction"]),
        "create_specimen_measure": lambda api: api.create_specimen_measure("2600000001", measure),
        "create_specimen_measures": lambda api: api.create_specimen_measures("2600000001", [measure]),
        "link_sequencing_sample_extraction": lambda api: api.link_sequencing_sample_extraction(
            vg_objects["sequencing_sample_lookup_1"], "2600000001C"),
        "upload_file_file_type": lambda api: api.upload_file(vcf, path=None, file_type=CVO),
        "poll_upload_status": lambda api: api.poll_upload_status(uploaded_file_id=1),
        "wait_for_annotation": lambda api: api.wait_for_annotation(uploaded_file_id=1, sleep=lambda _: None),
        "download_annotated": lambda api: api.download_annotated(uploaded_file_id=1, dest_path=vcf + ".out"),
        "annotate_vcf": lambda api: api.annotate_vcf(vcf, sleep=lambda _: None),  # checked before uploading
    }


GATED = ["create_patient", "create_specimen", "create_extraction", "create_specimen_measure",
         "create_specimen_measures", "link_sequencing_sample_extraction", "upload_file_file_type",
         "poll_upload_status", "wait_for_annotation", "download_annotated", "annotate_vcf"]


@pytest.mark.parametrize("method", GATED)
@responses.activate
def test_gated_method_raises_on_legacy_under_error(api, capabilities_url, vg_objects, vcf, method):
    responses.add(responses.GET, capabilities_url, status=404)

    with pytest.raises(UnsupportedFeatureError) as exc_info:
        _gated_calls(vg_objects, vcf)[method](api)
    assert exc_info.value.capabilities is ServerCapabilities.LEGACY
    assert len(responses.calls) == 1  # only the probe - nothing posted


@pytest.mark.parametrize("method", GATED)
@responses.activate
def test_gated_method_skips_on_legacy_under_skip(skip_api, capabilities_url, vg_objects, vcf, method, caplog):
    responses.add(responses.GET, capabilities_url, status=404)

    with caplog.at_level(logging.WARNING):
        assert _gated_calls(vg_objects, vcf)[method](skip_api) is None
    assert "Skipping" in caplog.text
    assert len(responses.calls) == 1  # only the probe - nothing posted


@responses.activate
def test_upload_file_metadata_raises_on_legacy_under_error(api, capabilities_url, vcf):
    responses.add(responses.GET, capabilities_url, status=404)

    with pytest.raises(UnsupportedFeatureError):
        api.upload_file(vcf, path=None, metadata={"genome_build": "GRCh37"})
    assert len(responses.calls) == 1  # only the probe - nothing posted


@responses.activate
def test_upload_file_metadata_dropped_on_legacy_under_skip(skip_api, server, capabilities_url, vcf, caplog):
    """ The file still goes up, as it did before metadata existed - only the metadata is skipped """
    responses.add(responses.GET, capabilities_url, status=404)
    upload_url = f"{server}/upload/api/v1/file_upload"
    responses.add(responses.POST, upload_url, json={"uploaded_file_id": 1}, status=200)

    with caplog.at_level(logging.WARNING):
        result = skip_api.upload_file(vcf, path=None, metadata={"genome_build": "GRCh37", "extraction": "2600000001C"})
    assert result == {"uploaded_file_id": 1}
    assert "Skipping: upload metadata" in caplog.text
    upload_call = responses.calls[-1]
    assert upload_call.request.url == upload_url  # no query params at all


@responses.activate
def test_gated_method_posts_when_supported(api, server, capabilities_url, capabilities_json, vg_objects):
    responses.add(responses.GET, capabilities_url, json=capabilities_json, status=200)
    responses.add(responses.POST, f"{server}/patients/api/v1/patient/", json={"id": 1}, status=201)

    assert api.create_patient(vg_objects["patient"]) == {"id": 1}


# ------------------------------------------------------------------ #
# upload_file(file_type=...)                                           #
# ------------------------------------------------------------------ #

@responses.activate
def test_upload_file_type_posts_when_accepted(skip_api, server, capabilities_url, capabilities_json, vcf):
    responses.add(responses.GET, capabilities_url, json=capabilities_json, status=200)
    upload_url = f"{server}/upload/api/v1/file_upload"
    responses.add(responses.POST, upload_url, json={"uploaded_file_id": 1}, status=200)

    assert skip_api.upload_file(vcf, path=None, file_type=CVO) == {"uploaded_file_id": 1}
    upload_call = responses.calls[-1]
    assert upload_call.request.url.startswith(upload_url)
    assert "file_type" not in upload_call.request.url  # a client-side gate only, not sent


@responses.activate
def test_upload_file_type_skipped_when_server_has_endpoint_but_not_type(skip_api, capabilities_url, capabilities_json,
                                                                        vcf):
    """ A VG3 server with the capabilities endpoint reports its file types but no CVO importer """
    responses.add(responses.GET, capabilities_url, status=200,
                  json={"version": "3.x", "git_hash": None, "features": [], "upload_file_types": ["vcf", "gene_list"]})

    assert skip_api.upload_file(vcf, path=None, file_type=CVO) is None
    assert len(responses.calls) == 1
