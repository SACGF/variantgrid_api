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
import warnings
from typing import Any, List, Optional, Tuple


class MockVariantGridAPI:
    """Drop-in replacement for VariantGridAPI that records all calls."""

    def __init__(self):
        # List of (method_name, args, kwargs) in call order
        self.calls: List[Tuple[str, tuple, dict]] = []
        self._return_values: dict = {}

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

    def upload_file(self, filename):
        self._record("upload_file", filename)
        return self._ret("upload_file", {"path": filename, "status": "ok"})

    def sequencing_run_has_vcf(self, sequencing_run, path=None):
        self._record("sequencing_run_has_vcf", sequencing_run, path)
        return self._ret("sequencing_run_has_vcf", True)

    def sequencing_run_name_has_vcf(self, sequencing_run_name, path=None):
        self._record("sequencing_run_name_has_vcf", sequencing_run_name, path)
        return self._ret("sequencing_run_name_has_vcf", True)
