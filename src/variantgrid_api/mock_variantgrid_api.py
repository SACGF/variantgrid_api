"""
Mock replacement for variantgrid_api.api_client.VariantGridAPI.

Records every call so tests can assert on what was sent without making any

Usage::

    api = MockVariantGridAPI()
    api.create_experiment("HAEM_25_080")
    api.assert_called_once("create_experiment")
    api.assert_not_called("sequencing_run_has_vcf")
"""
from __future__ import annotations

import copy
import logging
import warnings
from pathlib import Path
from typing import Any, List, Optional, Tuple, Union

from variantgrid_api.api_client import UnsupportedFeaturePolicy, UnsupportedFeatureError
from variantgrid_api.data_models import reference_json, ServerCapabilities, ServerFeature, UploadFileType

_UNSET = object()

# Every feature and upload file type the real client gates on - the shape of a current (VG4) server
MOCK_CAPABILITIES = ServerCapabilities(
    version="mock",
    features=frozenset(ServerFeature),
    upload_file_types=frozenset({UploadFileType.VCF, UploadFileType.GENE_COVERAGE,
                                 UploadFileType.DRAGEN_TSO500_ALL_FUSIONS,
                                 UploadFileType.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT,
                                 UploadFileType.GENE_LEVEL_CNV_VCF}),
)


class MockVariantGridAPI:
    """Drop-in replacement for VariantGridAPI that records all calls.

    capabilities defaults to MOCK_CAPABILITIES (a current server). Pass ServerCapabilities.LEGACY for an
    older server: gated calls then follow unsupported_feature_policy exactly as the real client does, and
    a skipped call is not recorded."""

    def __init__(self, capabilities: Optional[ServerCapabilities] = None,
                 unsupported_feature_policy=UnsupportedFeaturePolicy.ERROR,
                 logger: Optional[logging.Logger] = None):
        # List of (method_name, args, kwargs) in call order
        self.calls: List[Tuple[str, tuple, dict]] = []
        self._return_values: dict = {}
        self.capabilities = capabilities if capabilities is not None else MOCK_CAPABILITIES
        self.unsupported_feature_policy = unsupported_feature_policy
        self.logger = logger or logging.getLogger(__name__)

    # ------------------------------------------------------------------ #
    # Introspection helpers                                               #
    # ------------------------------------------------------------------ #

    def _record(self, method: str, *args, **kwargs):
        self.calls.append((method, args, kwargs))

    def _ret(self, method: str, default: Any) -> Any:
        return copy.deepcopy(self._return_values.get(method, default))

    def set_return(self, method: str, value: Any) -> None:
        """Override the default return value for a given method."""
        self._return_values[method] = value

    def get_calls(self, method: str) -> List[Tuple[tuple, dict]]:
        """Return all (args, kwargs) pairs recorded for *method*."""
        return [(a, kw) for name, a, kw in self.calls if name == method]

    def call_count(self, method: str) -> int:
        return len(self.get_calls(method))

    def assert_called_once(self, method: str) -> None:
        n = self.call_count(method)
        assert n == 1, f"Expected {method!r} called once, was called {n} times"

    def assert_not_called(self, method: str) -> None:
        n = self.call_count(method)
        assert n == 0, f"Expected {method!r} not called, was called {n} times"

    def reset(self) -> None:
        self.calls.clear()

    # ------------------------------------------------------------------ #
    # Capabilities — same gating as VariantGridAPI                        #
    # ------------------------------------------------------------------ #

    def supports(self, feature: Union[ServerFeature, str]) -> bool:
        return feature in self.capabilities.features

    def accepts_upload(self, file_type: Union[UploadFileType, str]) -> bool:
        return file_type in self.capabilities.upload_file_types

    def _unsupported(self, message: str) -> bool:
        message = f"{message} (server version '{self.capabilities.version}')"
        if self.unsupported_feature_policy == UnsupportedFeaturePolicy.SKIP:
            self.logger.warning("Skipping: %s", message)
            return False
        raise UnsupportedFeatureError(message, self.capabilities)

    def _require(self, feature: Union[ServerFeature, str]) -> bool:
        return self.supports(feature) or self._unsupported(f"server doesn't support feature '{feature}'")

    def _require_upload(self, file_type: Union[UploadFileType, str]) -> bool:
        return self.accepts_upload(file_type) or self._unsupported(f"server doesn't accept upload file type '{file_type}'")

    # ------------------------------------------------------------------ #
    # API surface — mirrors VariantGridAPI public methods exactly         #
    # ------------------------------------------------------------------ #

    def create_experiment(self, experiment):
        self._record("create_experiment", experiment)
        return self._ret("create_experiment", {"name": experiment})

    def create_enrichment_kit(self, enrichment_kit):
        self._record("create_enrichment_kit", enrichment_kit)
        return self._ret("create_enrichment_kit", {
            "name": enrichment_kit.name,
            "version": enrichment_kit.version,
            "gene_list": {
                "genelistgenesymbol_set": [
                    {"gene_symbol": "BRCA1"},
                    {"gene_symbol": "BRCA2"},
                ]
            },
        })

    def create_sequencer_model(self, sequencer_model):
        self._record("create_sequencer_model", sequencer_model)
        return self._ret("create_sequencer_model", {"model": sequencer_model.model})

    def create_sequencer(self, sequencer):
        self._record("create_sequencer", sequencer)
        return self._ret("create_sequencer", {"name": sequencer.name})

    def create_sequencing_run(self, sequencing_run):
        self._record("create_sequencing_run", sequencing_run)
        return self._ret("create_sequencing_run", {"name": sequencing_run.name})

    def create_sample_sheet(self, sample_sheet):
        self._record("create_sample_sheet", sample_sheet)
        return self._ret("create_sample_sheet", {"path": sample_sheet.path})

    def create_joint_called_vcf(self, joint_called_vcf):
        self._record("create_joint_called_vcf", joint_called_vcf)
        return self._ret("create_joint_called_vcf", {"path": joint_called_vcf.path})

    def create_sample_sheet_combined_vcf_file(self, sscvf):
        """Deprecated alias for :meth:`create_joint_called_vcf` - use that instead. """
        warnings.warn(
            "create_sample_sheet_combined_vcf_file is deprecated; use create_joint_called_vcf instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        self._record("create_sample_sheet_combined_vcf_file", sscvf)
        return self._ret("create_sample_sheet_combined_vcf_file", {"path": sscvf.path})

    def create_sequencing_data(self, sample_sheet_lookup, sequencing_files):
        self._record("create_sequencing_data", sample_sheet_lookup, sequencing_files)
        return self._ret("create_sequencing_data", {"created": len(sequencing_files)})

    def create_qc_gene_list(self, qc_gene_list):
        self._record("create_qc_gene_list", qc_gene_list)
        return self._ret("create_qc_gene_list", {})

    def create_multiple_qc_gene_lists(self, qc_gene_lists):
        self._record("create_multiple_qc_gene_lists", qc_gene_lists)
        return self._ret("create_multiple_qc_gene_lists", {"created": len(qc_gene_lists)})

    def create_qc_exec_stats(self, qc_exec_stats):
        self._record("create_qc_exec_stats", qc_exec_stats)
        return self._ret("create_qc_exec_stats", {})

    def create_multiple_qc_exec_stats(self, qc_exec_stats):
        self._record("create_multiple_qc_exec_stats", qc_exec_stats)
        return self._ret("create_multiple_qc_exec_stats", {"created": len(qc_exec_stats)})

    def create_multiple_qc_gene_coverage(self, qc_gene_coverage_list):
        self._record("create_multiple_qc_gene_coverage", qc_gene_coverage_list)
        return self._ret("create_multiple_qc_gene_coverage", {"created": len(qc_gene_coverage_list)})

    def create_patient(self, patient):
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._record("create_patient", patient)
        return self._ret("create_patient", {"id": 1, **patient.to_dict()})

    def create_specimen(self, specimen):
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._record("create_specimen", specimen)
        return self._ret("create_specimen", {"id": 1, **specimen.to_dict()})

    def create_extraction(self, extraction):
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._record("create_extraction", extraction)
        return self._ret("create_extraction", {"id": 1, **extraction.to_dict()})

    def create_specimen_measure(self, specimen_reference, measure):
        if not self._require(ServerFeature.SPECIMEN_MEASURES):
            return None
        self._record("create_specimen_measure", specimen_reference, measure)
        return self._ret("create_specimen_measure", {"id": 1, "specimen": reference_json(specimen_reference),
                                                     **measure.to_dict()})

    def create_specimen_measures(self, specimen_reference, measures):
        if not self._require(ServerFeature.SPECIMEN_MEASURES):
            return None
        self._record("create_specimen_measures", specimen_reference, measures)
        return self._ret("create_specimen_measures", {"specimen": reference_json(specimen_reference),
                                                      "measures": [m.to_dict() for m in measures]})

    def link_sequencing_sample_extraction(self, sequencing_sample_lookup, extraction_reference):
        if not self._require(ServerFeature.LINK_EXTRACTION):
            return None
        self._record("link_sequencing_sample_extraction", sequencing_sample_lookup, extraction_reference)
        return self._ret("link_sequencing_sample_extraction", {
            "sequencing_sample": sequencing_sample_lookup.sample_name,
            "match_status": "Matched",
            "match_error": None,
            "extraction": str(reference_json(extraction_reference)),
        })

    def upload_file(self, filename, path=_UNSET, metadata=None, file_type=None):
        if metadata and not self.supports(ServerFeature.UPLOAD_METADATA):
            self._unsupported(f"upload metadata for '{filename}' (server doesn't support feature 'upload_metadata')")
            metadata = None
        if file_type and not self._require_upload(file_type):
            return None
        if path is _UNSET:
            path = filename
        # metadata / file_type only recorded when sent, so existing assertions on the recorded kwargs still hold
        extra = {"metadata": metadata} if metadata is not None else {}
        if file_type is not None:
            extra["file_type"] = file_type
        self._record("upload_file", filename, path=path, **extra)
        return self._ret("upload_file", {"uploaded_file_id": 1, "sha256_hash": "deadbeef",
                                         "path": path, "status": "ok"})

    def poll_upload_status(self, uploaded_file_id=None, sha256=None):
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        self._record("poll_upload_status", uploaded_file_id, sha256)
        return self._ret("poll_upload_status", {
            "uploaded_file_id": uploaded_file_id,
            "sha256_hash": sha256,
            "annotation_complete": True,
            "error": None,
        })

    def wait_for_annotation(self, uploaded_file_id=None, sha256=None,
                            timeout=3600, poll_interval=10, sleep=None, max_transient_errors=5):
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        self._record("wait_for_annotation", uploaded_file_id, sha256,
                     timeout=timeout, poll_interval=poll_interval, sleep=sleep,
                     max_transient_errors=max_transient_errors)
        return self._ret("wait_for_annotation", {
            "uploaded_file_id": uploaded_file_id,
            "annotation_complete": True,
            "error": None,
        })

    def download_annotated(self, uploaded_file_id=None, sha256=None, export_type="vcf",
                           dest_path=None, timeout=3600, poll_interval=10, sleep=None):
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        self._record("download_annotated", uploaded_file_id, sha256, export_type=export_type,
                     dest_path=dest_path, timeout=timeout, poll_interval=poll_interval, sleep=sleep)
        default = Path(dest_path) if dest_path is not None else Path(f"download.{export_type}")
        return self._ret("download_annotated", default)

    def annotate_vcf(self, filename, export_type="vcf", dest_path=None,
                     timeout=3600, poll_interval=10, sleep=None):
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        self._record("annotate_vcf", filename, export_type=export_type, dest_path=dest_path,
                     timeout=timeout, poll_interval=poll_interval, sleep=sleep)
        default = Path(dest_path) if dest_path is not None else Path(f"download.{export_type}")
        return self._ret("annotate_vcf", default)

    def sequencing_run_has_vcf(self, sequencing_run, path=None):
        self._record("sequencing_run_has_vcf", sequencing_run, path)
        return self._ret("sequencing_run_has_vcf", True)

    def sequencing_run_name_has_vcf(self, sequencing_run_name, path=None):
        self._record("sequencing_run_name_has_vcf", sequencing_run_name, path)
        return self._ret("sequencing_run_name_has_vcf", True)
