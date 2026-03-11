import pytest
from variantgrid_api.data_models import SequencerModel


@pytest.mark.parametrize("name,expected_model,expected_convention", [
    ("M02027",  "MiSeq",        "M"),
    ("H001",    "HiSeq",        "H"),
    ("NB551234","NextSeq",      "M"),
    ("NS500123","NextSeq",      "M"),
    ("A01234",  "NovaSeq 6000", "M"),
    ("LH00999", "NovaSeq X",    "M"),
    ("XYZ999",  "Unknown",      "U"),
])
def test_from_sequencer_name(name, expected_model, expected_convention):
    sm = SequencerModel.from_sequencer_name(name)
    assert sm.model == expected_model
    assert sm.data_naming_convention == expected_convention
    assert sm.manufacturer.name == "Illumina"


def test_from_sequencer_name_matches_conftest_example():
    # The conftest uses sequencer "SN1101" which doesn't match any prefix → Unknown
    # But the run name contains "M02027" as the instrument — verify MiSeq detection
    sm = SequencerModel.from_sequencer_name("M02027")
    assert sm.model == "MiSeq"


def test_lh_prefix_not_matched_as_h():
    """NovaSeq X instruments start with LH — must not be classified as HiSeq."""
    sm = SequencerModel.from_sequencer_name("LH00187")
    assert sm.model == "NovaSeq X"
    assert sm.data_naming_convention == "M"
