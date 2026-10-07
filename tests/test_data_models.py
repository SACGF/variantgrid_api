import pytest

from variantgrid_api.data_models import (
    SequencingRun, SampleSheetLookup, JointCalledVCF, SampleSheetCombinedVCFFile,
    SingleSampleVCF, VCFFile, VariantCaller,
)

def test_get_date_from_name(vg_objects):
    d = SequencingRun.get_date_from_name(vg_objects["SEQUENCING_RUN_NAME"])
    assert d is not None

def test_sample_sheet_lookup_from_sample_sheet(vg_objects):
    ss = vg_objects["sample_sheet"]
    ssl = SampleSheetLookup.from_sample_sheet(ss)
    d = ssl.to_dict()
    assert d


def test_sample_sheet_combined_vcf_file_is_deprecated_alias(vg_objects):
    """SampleSheetCombinedVCFFile still works, is a JointCalledVCF, but warns on use."""
    lookup = vg_objects["sample_sheet_lookup"]
    variant_caller = VariantCaller(name="VarDict", version="1.8.2")
    with pytest.warns(DeprecationWarning):
        obj = SampleSheetCombinedVCFFile(
            path="/data/combined.vcf.gz", sample_sheet_lookup=lookup, variant_caller=variant_caller
        )
    assert isinstance(obj, JointCalledVCF)
    # Wire format unchanged - serialises the same as the canonical class
    canonical = JointCalledVCF(
        path="/data/combined.vcf.gz", sample_sheet_lookup=lookup, variant_caller=variant_caller
    )
    assert obj.to_dict() == canonical.to_dict()


def test_vcf_file_is_deprecated_alias():
    """VCFFile still works, is a SingleSampleVCF, but warns on use."""
    variant_caller = VariantCaller(name="GATK", version="4.1.9.0")
    with pytest.warns(DeprecationWarning):
        obj = VCFFile(path="/data/sample.vcf.gz", variant_caller=variant_caller)
    assert isinstance(obj, SingleSampleVCF)
    canonical = SingleSampleVCF(path="/data/sample.vcf.gz", variant_caller=variant_caller)
    assert obj.to_dict() == canonical.to_dict()
