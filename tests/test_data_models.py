from variantgrid_api.data_models import SequencingRun, SampleSheetLookup

def test_get_date_from_name(vg_objects):
    d = SequencingRun.get_date_from_name(vg_objects["SEQUENCING_RUN_NAME"])
    assert d is not None

def test_sample_sheet_lookup_from_sample_sheet(vg_objects):
    ss = vg_objects["sample_sheet"]
    ssl = SampleSheetLookup.from_sample_sheet(ss)
    d = ssl.to_dict()
    assert d
