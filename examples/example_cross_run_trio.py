"""
    A family trio joint-called from samples that were sequenced on two different runs.

    Each run's sample sheet goes up as usual, then one JointCalledVCF names its three members
    explicitly via SequencingSampleLookup. The joint call's own sample_sheet_lookup stays the
    owning run - the one the VCF path sits under.
"""
import argparse
import os

from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit, SequencerModel, Sequencer, SequencingRun, SequencingSample, \
    SampleSheet, JointCalledVCF, VariantCaller, SampleSheetLookup, SequencingSampleLookup, Manufacturer


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client - cross-run family trio")
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--api-token', required=True, help='API token for authentication')
    parser.add_argument('--step', required=False, help='Run a single step (default: run all)')
    return parser.parse_args()


def _make_run_and_sample_sheet(data_dir, run_name, sample_names, sheet_hash, enrichment_kit, sequencer):
    seq_run_dir = os.path.join(data_dir, "wgs", run_name)
    sequencing_run = SequencingRun(path=seq_run_dir,
                                   name=run_name,
                                   date=SequencingRun.get_date_from_name(run_name),
                                   sequencer=sequencer.name,
                                   experiment="WGS_26_033",
                                   enrichment_kit=enrichment_kit)

    sequencing_samples = [
        SequencingSample(sample_id=sample_name,
                         sample_number=i,
                         barcode="GCCAAT",
                         enrichment_kit=enrichment_kit)
        for i, sample_name in enumerate(sample_names, start=1)
    ]

    sample_sheet = SampleSheet(path=os.path.join(seq_run_dir, "SampleSheet.csv"),
                               sequencing_run=sequencing_run,
                               file_last_modified=1725941707.0033002,
                               hash=sheet_hash,
                               sequencing_samples=sequencing_samples)
    return sequencing_run, sample_sheet


def test_api(server, api_token, step=None):
    proj_dir = os.path.dirname(os.path.dirname(__file__))  # In "examples"
    data_dir = os.path.join(proj_dir, "tests", 'test_data')

    enrichment_kit = EnrichmentKit(name='wgs', version=1)
    manufacturer = Manufacturer(name="Illumina")
    sequencer_model = SequencerModel(model="NovaSeq X", manufacturer=manufacturer, data_naming_convention="M")
    sequencer = Sequencer(name="LH00639", sequencer_model=sequencer_model)

    # The proband was sequenced on one run, the parents on another
    proband_run, proband_sheet = _make_run_and_sample_sheet(
        data_dir, "WGS_26_033_260508_LH00639_0099_A23LFGJLT3", ["PND151_AFF"],
        "f0ac87bcae3f0e56b3f65b70fd6389ce", enrichment_kit, sequencer)
    parents_run, parents_sheet = _make_run_and_sample_sheet(
        data_dir, "WGS_26_021_260327_LH00639_0087_A22KDEFGH1", ["PND151_MUM", "PND151_DAD"],
        "b31c0e5a6dbb8f9c1e4f7a2d3c5b6e80", enrichment_kit, sequencer)

    proband_sheet_lookup = SampleSheetLookup.from_sample_sheet(proband_sheet)
    parents_sheet_lookup = SampleSheetLookup.from_sample_sheet(parents_sheet)

    # The trio VCF lives under the proband's run, so that is the owning sample sheet
    trio_vcf_filename = os.path.join(proband_run.path, "CLIN_WGS_PND151_family", "VCFs",
                                     "CLIN_WGS_PND151.vcf.gz")
    trio_vcf = JointCalledVCF(
        path=trio_vcf_filename,
        sample_sheet_lookup=proband_sheet_lookup,
        variant_caller=VariantCaller(name="GATK", version="4.1.9.0"),
        sequencing_samples=[
            SequencingSampleLookup(sample_sheet_lookup=proband_sheet_lookup, sample_name="PND151_AFF"),
            SequencingSampleLookup(sample_sheet_lookup=parents_sheet_lookup, sample_name="PND151_MUM"),
            SequencingSampleLookup(sample_sheet_lookup=parents_sheet_lookup, sample_name="PND151_DAD"),
        ])

    #########################
    # Call API

    vg_api = VariantGridAPI(server, api_token)

    # Both sample sheets go up before the trio VCF, so every member lookup resolves
    API_STEPS = {
        "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
        "sequencer_model": lambda: vg_api.create_sequencer_model(sequencer_model),
        "sequencer": lambda: vg_api.create_sequencer(sequencer),
        "proband_sequencing_run": lambda: vg_api.create_sequencing_run(proband_run),
        "proband_sample_sheet": lambda: vg_api.create_sample_sheet(proband_sheet),
        "parents_sequencing_run": lambda: vg_api.create_sequencing_run(parents_run),
        "parents_sample_sheet": lambda: vg_api.create_sample_sheet(parents_sheet),
        "trio_joint_called_vcf": lambda: vg_api.create_joint_called_vcf(trio_vcf),
        "upload_trio_vcf_file": lambda: vg_api.upload_file(trio_vcf_filename),
    }

    for name, func in API_STEPS.items():
        if step:
            if name != step:
                continue
        print(f"{name=}")
        result = func()
        print(f"{result=}")
        print("-" * 50)


if __name__ == "__main__":
    args = parse_args()
    test_api(args.server, args.api_token,
             step=args.step)
