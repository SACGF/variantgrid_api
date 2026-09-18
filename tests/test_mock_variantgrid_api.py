"""Tests for MockVariantGridAPI — verifies call recording, return-value
overrides, and assertion helpers work correctly."""
import pytest

from variantgrid_api.api_client import UnsupportedFeaturePolicy, UnsupportedFeatureError
from variantgrid_api.data_models import ServerCapabilities, ServerFeature, UploadFileType
from variantgrid_api.mock_variantgrid_api import MockVariantGridAPI


@pytest.fixture
def mock_api():
    return MockVariantGridAPI()


# ------------------------------------------------------------------ #
# Call recording                                                       #
# ------------------------------------------------------------------ #

def test_calls_are_recorded(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    assert mock_api.call_count("create_experiment") == 1
    assert mock_api.call_count("create_sequencing_run") == 0


def test_multiple_calls_accumulate(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.create_experiment(vg_objects["experiment"])
    assert mock_api.call_count("create_experiment") == 2


def test_get_calls_returns_args(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    calls = mock_api.get_calls("create_experiment")
    assert len(calls) == 1
    args, kwargs = calls[0]
    assert args[0] == vg_objects["experiment"]


def test_calls_from_different_methods_are_independent(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.create_sequencing_run(vg_objects["sequencing_run"])
    assert mock_api.call_count("create_experiment") == 1
    assert mock_api.call_count("create_sequencing_run") == 1
    assert mock_api.call_count("create_sample_sheet") == 0


# ------------------------------------------------------------------ #
# Assertion helpers                                                    #
# ------------------------------------------------------------------ #

def test_assert_called_once_passes(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.assert_called_once("create_experiment")  # should not raise


def test_assert_called_once_fails_when_not_called(mock_api):
    with pytest.raises(AssertionError):
        mock_api.assert_called_once("create_experiment")


def test_assert_called_once_fails_when_called_twice(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.create_experiment(vg_objects["experiment"])
    with pytest.raises(AssertionError):
        mock_api.assert_called_once("create_experiment")


def test_assert_not_called_passes(mock_api):
    mock_api.assert_not_called("create_experiment")  # should not raise


def test_assert_not_called_fails_when_called(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    with pytest.raises(AssertionError):
        mock_api.assert_not_called("create_experiment")


# ------------------------------------------------------------------ #
# Default return values                                                #
# ------------------------------------------------------------------ #

def test_default_return_create_experiment(mock_api, vg_objects):
    result = mock_api.create_experiment(vg_objects["experiment"])
    assert result == {"name": vg_objects["experiment"]}


def test_default_return_sequencing_run_has_vcf(mock_api, vg_objects):
    result = mock_api.sequencing_run_has_vcf(vg_objects["sequencing_run"])
    assert result is True


def test_default_return_create_sequencing_data(mock_api, vg_objects):
    result = mock_api.create_sequencing_data(
        vg_objects["sample_sheet_lookup"], vg_objects["sequencing_files"]
    )
    assert result == {"created": len(vg_objects["sequencing_files"])}


def test_default_return_create_multiple_qc_gene_lists(mock_api, vg_objects):
    result = mock_api.create_multiple_qc_gene_lists(vg_objects["qc_gene_lists"])
    assert result == {"created": len(vg_objects["qc_gene_lists"])}


def test_create_joint_called_vcf_records_and_returns(mock_api, vg_objects):
    jcv = vg_objects["joint_called_vcf"]
    result = mock_api.create_joint_called_vcf(jcv)
    mock_api.assert_called_once("create_joint_called_vcf")
    assert result == {"path": jcv.path}


def test_create_sample_sheet_combined_vcf_file_deprecated_alias(mock_api, vg_objects):
    jcv = vg_objects["joint_called_vcf"]
    with pytest.warns(DeprecationWarning):
        result = mock_api.create_sample_sheet_combined_vcf_file(jcv)
    mock_api.assert_called_once("create_sample_sheet_combined_vcf_file")
    assert result == {"path": jcv.path}


# ------------------------------------------------------------------ #
# Uploaded file annotation flow                                        #
# ------------------------------------------------------------------ #

def test_mock_poll_upload_status_records_and_returns(mock_api):
    result = mock_api.poll_upload_status(uploaded_file_id=123)
    mock_api.assert_called_once("poll_upload_status")
    assert result["annotation_complete"] is True


def test_mock_wait_for_annotation_records_and_returns(mock_api):
    result = mock_api.wait_for_annotation(uploaded_file_id=123)
    mock_api.assert_called_once("wait_for_annotation")
    assert result["annotation_complete"] is True


def test_mock_download_annotated_defaults_to_dest_path(mock_api):
    result = mock_api.download_annotated(uploaded_file_id=123, dest_path="/tmp/out.vcf.gz")
    mock_api.assert_called_once("download_annotated")
    assert str(result) == "/tmp/out.vcf.gz"


def test_mock_annotate_vcf_records_and_returns(mock_api):
    result = mock_api.annotate_vcf("input.vcf", export_type="csv")
    mock_api.assert_called_once("annotate_vcf")
    assert str(result) == "download.csv"


# ------------------------------------------------------------------ #
# set_return overrides                                                 #
# ------------------------------------------------------------------ #

def test_set_return_overrides_default(mock_api, vg_objects):
    mock_api.set_return("create_experiment", {"name": "overridden", "id": 99})
    result = mock_api.create_experiment(vg_objects["experiment"])
    assert result == {"name": "overridden", "id": 99}


def test_set_return_is_deep_copied(mock_api, vg_objects):
    sentinel = {"name": "exp", "extra": []}
    mock_api.set_return("create_experiment", sentinel)
    r1 = mock_api.create_experiment(vg_objects["experiment"])
    r1["extra"].append("mutated")
    r2 = mock_api.create_experiment(vg_objects["experiment"])
    assert r2["extra"] == [], "set_return should return independent copies each call"


def test_set_return_sequencing_run_has_vcf_false(mock_api, vg_objects):
    mock_api.set_return("sequencing_run_has_vcf", False)
    assert mock_api.sequencing_run_has_vcf(vg_objects["sequencing_run"]) is False


# ------------------------------------------------------------------ #
# reset                                                                #
# ------------------------------------------------------------------ #

def test_reset_clears_calls(mock_api, vg_objects):
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.reset()
    assert mock_api.call_count("create_experiment") == 0


def test_reset_preserves_return_overrides(mock_api, vg_objects):
    mock_api.set_return("create_experiment", {"id": 7})
    mock_api.create_experiment(vg_objects["experiment"])
    mock_api.reset()
    result = mock_api.create_experiment(vg_objects["experiment"])
    assert result == {"id": 7}


# ------------------------------------------------------------------ #
# Patient / specimen / extraction, measures, links                     #
# ------------------------------------------------------------------ #

def test_mock_create_patient_specimen_extraction(mock_api, vg_objects):
    assert mock_api.create_patient(vg_objects["patient"])["patient_code"] == "C0000001"
    assert mock_api.create_specimen(vg_objects["specimen"])["reference_id"] == "2600000001"
    assert mock_api.create_extraction(vg_objects["dna_extraction"])["reference_id"] == "2600000001C"
    for method in ("create_patient", "create_specimen", "create_extraction"):
        mock_api.assert_called_once(method)


def test_mock_create_specimen_measure_records_args(mock_api, vg_objects):
    measure = vg_objects["specimen_measures"][0]
    result = mock_api.create_specimen_measure("2600000001", measure)
    assert mock_api.get_calls("create_specimen_measure") == [(("2600000001", measure), {})]
    assert result["specimen"] == "2600000001"


def test_mock_create_specimen_measures(mock_api, vg_objects):
    result = mock_api.create_specimen_measures("2600000001", vg_objects["specimen_measures"])
    mock_api.assert_called_once("create_specimen_measures")
    assert len(result["measures"]) == len(vg_objects["specimen_measures"])


def test_mock_link_sequencing_sample_extraction_pending_override(mock_api, vg_objects):
    mock_api.set_return("link_sequencing_sample_extraction", {"match_status": "Pending"})
    result = mock_api.link_sequencing_sample_extraction(vg_objects["sequencing_sample_lookup_1"], "2600000001C")
    assert result == {"match_status": "Pending"}
    mock_api.assert_called_once("link_sequencing_sample_extraction")


def test_mock_upload_file_records_metadata_only_when_sent(mock_api):
    mock_api.upload_file("a.vcf")
    mock_api.upload_file("b.vcf", path=None, metadata={"genome_build": "GRCh37"})
    calls = mock_api.get_calls("upload_file")
    assert calls[0] == (("a.vcf",), {"path": "a.vcf"})
    assert calls[1] == (("b.vcf",), {"path": None, "metadata": {"genome_build": "GRCh37"}})


# ------------------------------------------------------------------ #
# Capabilities                                                         #
# ------------------------------------------------------------------ #

def test_mock_default_capabilities_support_everything_gated(mock_api, vg_objects):
    assert mock_api.supports("patients")
    assert mock_api.accepts_upload("dragen_tso500_combined_variant_output")
    mock_api.create_patient(vg_objects["patient"])
    mock_api.upload_file("cvo.tsv", path=None, file_type="dragen_tso500_combined_variant_output")
    mock_api.assert_called_once("create_patient")
    assert mock_api.get_calls("upload_file") == [
        (("cvo.tsv",), {"path": None, "file_type": "dragen_tso500_combined_variant_output"})]


def test_mock_default_capabilities_accept_enums(mock_api):
    assert all(mock_api.supports(feature) for feature in ServerFeature)
    assert mock_api.supports("upload_status")
    assert mock_api.accepts_upload(UploadFileType.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT)
    assert not mock_api.accepts_upload(UploadFileType.PED)
    mock_api.upload_file("cvo.tsv", path=None, file_type=UploadFileType.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT)
    mock_api.assert_called_once("upload_file")


def test_mock_legacy_skips_under_skip(vg_objects):
    mock_api = MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY,
                                  unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)
    assert not mock_api.supports("patients")
    assert mock_api.create_patient(vg_objects["patient"]) is None
    assert mock_api.link_sequencing_sample_extraction(vg_objects["sequencing_sample_lookup_1"], "2600000001C") is None
    assert mock_api.upload_file("cvo.tsv", path=None, file_type="dragen_tso500_combined_variant_output") is None
    mock_api.assert_not_called("create_patient")
    mock_api.assert_not_called("link_sequencing_sample_extraction")
    mock_api.assert_not_called("upload_file")

    # Metadata is dropped but the file still uploads
    assert mock_api.upload_file("a.vcf", path=None, metadata={"genome_build": "GRCh37"}) is not None
    assert mock_api.get_calls("upload_file") == [(("a.vcf",), {"path": None})]


def test_mock_legacy_raises_under_error(vg_objects):
    mock_api = MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY)
    with pytest.raises(UnsupportedFeatureError):
        mock_api.create_specimen_measures("2600000001", vg_objects["specimen_measures"])
    mock_api.assert_not_called("create_specimen_measures")


def test_mock_legacy_skips_annotation_flow():
    mock_api = MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY,
                                  unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)
    assert mock_api.poll_upload_status(uploaded_file_id=1) is None
    assert mock_api.annotate_vcf("input.vcf") is None
    mock_api.assert_not_called("poll_upload_status")
    mock_api.assert_not_called("annotate_vcf")
