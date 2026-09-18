"""
    A BAM-first run - the sequencer emits BAM directly (or FastQs aren't kept), so there are no FastQs.

    Same shape as example_haem_20_999.py, but SequencingFile has no fastq_r1/fastq_r2. The record sent
    is just bam_file + vcf_file, and the server resolves the sample from sample_name + the sample sheet.
"""
import argparse
import dataclasses
import os
from typing import List, Dict

from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit, SequencerModel, Sequencer, SequencingRun, SequencingSample, \
    SampleSheet, JointCalledVCF, VariantCaller, SampleSheetLookup, Aligner, SingleSampleVCF, BamFile, SequencingFile, \
    SequencingSampleLookup, QC, QCGeneList, Manufacturer


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client - BAM-first run (no FastQs)")
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--api-token', required=True, help='API token for authentication')
    parser.add_argument('--step', required=False, help='Run a single step (default: run all)')
    return parser.parse_args()


def _get_qc_by_sample_name(sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]) -> Dict[str, QC]:
    qc_by_name = {}
    for sf in sequencing_files:
        sequencing_sample_lookup = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup,
                                                          sample_name=sf.sample_name)
        qc_by_name[sf.sample_name] = QC(sequencing_sample_lookup=sequencing_sample_lookup,
                                        bam_file=dataclasses.replace(sf.bam_file, aligner=None),
                                        vcf_file=dataclasses.replace(sf.get_vcf_files()[0], variant_caller=None))
    return qc_by_name


def test_api(server, api_token, step=None):
    proj_dir = os.path.dirname(os.path.dirname(__file__))  # In "examples"
    data_dir = os.path.join(proj_dir, "tests", 'test_data')
    SEQUENCING_RUN_NAME = "Haem_20_998_201230_M02027_0111_000000000_JFT78"
    seq_run_dir = os.path.join(data_dir, "idt_haem", SEQUENCING_RUN_NAME)

    def seq_run_path(path):
        return os.path.join(seq_run_dir, path)

    experiment = "HAEM_20_998"
    enrichment_kit = EnrichmentKit(name='idt_haem', version=1)

    manufacturer = Manufacturer(name="Illumina")
    sequencer_model = SequencerModel(model="HiSeq 2500", manufacturer=manufacturer, data_naming_convention="H")
    sequencer_name = "SN1101"
    sequencer = Sequencer(name=sequencer_name, sequencer_model=sequencer_model)
    sequencing_run = SequencingRun(path=seq_run_dir,
                                   name=SEQUENCING_RUN_NAME,
                                   date=SequencingRun.get_date_from_name(SEQUENCING_RUN_NAME),
                                   sequencer=sequencer_name,
                                   experiment=experiment,
                                   enrichment_kit=enrichment_kit)

    sample_names = ["fake_sample_1", "fake_sample_2"]
    sequencing_samples = [
        SequencingSample(sample_id=sample_name,
                         sample_number=i,
                         barcode="GCCAAT",
                         enrichment_kit=enrichment_kit)
        for i, sample_name in enumerate(sample_names, start=1)
    ]

    sample_sheet = SampleSheet(
        path=seq_run_path("SampleSheet.csv"),
        sequencing_run=sequencing_run,
        file_last_modified=1725941707.0033002,
        hash="6d1e3b0a9c2f4e5d8b7a6c5d4e3f2a1b",
        sequencing_samples=sequencing_samples)

    sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)

    joint_called_vcf = JointCalledVCF(
        path=seq_run_path(f"2_variants/{SEQUENCING_RUN_NAME}.vardict.hg38.vcf.gz"),
        sample_sheet_lookup=sample_sheet_lookup,
        variant_caller=VariantCaller(name="VarDict", version="1.8.2"))

    aligner = Aligner(name='BWA', version="0.7.18")
    variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")

    # No fastq_r1/fastq_r2 - only BAM + VCF
    sequencing_files = [
        SequencingFile(sample_name=sample_name,
                       bam_file=BamFile(path=seq_run_path(f"1_BAM/{sample_name}.hg38.bam"),
                                        aligner=aligner),
                       vcf_files=[SingleSampleVCF(path=seq_run_path(f"2_variants/gatk_per_sample/{sample_name}.gatk.hg38.vcf.gz"),
                                                  variant_caller=variant_caller_gatk)])
        for sample_name in sample_names
    ]

    # QC is matched against BAM/VCF, so works the same without FastQs
    qc_by_sample_name = _get_qc_by_sample_name(sample_sheet_lookup, sequencing_files)
    gene_list = ["TUBA1A", "TUBA8", "FLNA", "TUBB2B", "TUBB3", "COL4A1", "KIAA1279"]
    qc_gene_lists = [
        QCGeneList(path=seq_run_path(f"0_goi/{SEQUENCING_RUN_NAME}_{sample_name}.txt"),
                   qc=qc_by_sample_name[sample_name],
                   gene_list=gene_list)
        for sample_name in sample_names
    ]

    #########################
    # Call API

    vg_api = VariantGridAPI(server, api_token)

    API_STEPS = {
        "experiment": lambda: vg_api.create_experiment(experiment),
        "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
        "sequencer_model": lambda: vg_api.create_sequencer_model(sequencer_model),
        "sequencer": lambda: vg_api.create_sequencer(sequencer),
        "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
        "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
        "joint_called_vcf": lambda: vg_api.create_joint_called_vcf(joint_called_vcf),
        "sequencing_data": lambda: vg_api.create_sequencing_data(sample_sheet_lookup, sequencing_files),
        "qc_gene_lists": lambda: vg_api.create_multiple_qc_gene_lists(qc_gene_lists),
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
