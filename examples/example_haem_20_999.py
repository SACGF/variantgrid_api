import argparse
import dataclasses
import os
from datetime import datetime
from hashlib import md5
from typing import List, Dict

from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit, SequencingRun, SequencingSample, SampleSheet, \
    SampleSheetCombinedVCFFile, VariantCaller, SampleSheetLookup, Aligner, VCFFile, BamFile, SequencingFile, \
    SequencingSampleLookup, QC, QCGeneList, QCExecStats, QCGeneCoverage


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client")
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--api-token', required=True, help='API token for authentication')
    parser.add_argument('--step', required=False, help='Run a single step (default: run all)')
    return parser.parse_args()


def _file_md5sum(filename):
    m = md5()
    with open(filename, "rb") as f:
        m.update(f.read())
    return m.hexdigest()


def _get_qc_by_sample_name(sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]) -> Dict[str, QC]:
    bam_and_vcf_by_name = {}
    for sf in sequencing_files:
        bam_file = dataclasses.replace(sf.bam_file, aligner=None)
        vcf_file = dataclasses.replace(sf.vcf_file, variant_caller=None)
        bam_and_vcf_by_name[sf.sample_name] = (bam_file, vcf_file)

    qc_by_name = {}
    for sample_name, (bam_file, vcf_file) in bam_and_vcf_by_name.items():
        sequencing_sample_lookup = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup,
                                                          sample_name=sample_name)
        qc_by_name[sample_name] = QC(sequencing_sample_lookup=sequencing_sample_lookup,
                                     bam_file=bam_file,
                                     vcf_file=vcf_file)
    return qc_by_name



def test_api(server, api_token, step=None):
    data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
    SEQUENCING_RUN_NAME = "Haem_20_999_201231_M02027_0112_000000000_JFT79"
    seq_run_dir = os.path.join(data_dir, "idt_haem", SEQUENCING_RUN_NAME)

    def seq_run_path(path):
        return os.path.join(seq_run_dir, path)

    experiment = "HAEM_20_999"
    enrichment_kit = EnrichmentKit(name='idt_haem', version=1)
    seq_date = SequencingRun.get_date_from_name(SEQUENCING_RUN_NAME)
    sequencing_run = SequencingRun(path=seq_run_dir,
                                   name=SEQUENCING_RUN_NAME,
                                   date=seq_date,
                                   sequencer="SN1101",
                                   experiment=experiment,
                                   enrichment_kit=enrichment_kit)

    sequencing_samples = [
        SequencingSample(sample_id="fake_sample_1",
                         sample_project=None,
                         sample_number=1,
                         barcode="GCCAAT",
                         enrichment_kit=enrichment_kit,
                         is_control=False,
                         failed=False,
                         data=[
                             {
                                 "column": "SAPOrderNumber",
                                 "value": "SAP1000001"
                             }
                         ]),
        SequencingSample(sample_id="fake_sample_2",
                         sample_project=None,
                         sample_number=2,
                         barcode="CAGATC",
                         enrichment_kit=enrichment_kit,
                         is_control=False,
                         failed=False,
                         data=[
                             {
                                 "column": "SAPOrderNumber",
                                 "value": "SAP1000002"
                             }
                         ])
    ]

    sample_sheet = SampleSheet(
        path=seq_run_path("SampleSheet.csv"),
        sequencing_run=sequencing_run,
        file_last_modified=1725941707.0033002,
        hash="f0ac87bcae3f0e56b3f65b70fd6389ce",
        sequencing_samples=sequencing_samples)

    sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)

    variant_caller_var_dict = VariantCaller(name="VarDict", version="1.8.2")
    combo_vcf_filename = seq_run_path("2_variants/Haem_20_999_201231_M02027_0112_000000000_JFT79.vardict.hg38.vcf.gz")
    sample_sheet_combined_vcf_file = SampleSheetCombinedVCFFile(
        path=combo_vcf_filename,
        sample_sheet_lookup=sample_sheet_lookup,
        variant_caller=variant_caller_var_dict)

    aligner = Aligner(name='BWA', version="0.7.18")
    variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")
    bam_file_1 = BamFile(
        path=seq_run_path("1_BAM/fake_sample_1.hg38.bam"),
        aligner=aligner)
    vcf_file_1 = VCFFile(
        path=seq_run_path("2_variants/gatk_per_sample/fake_sample_1.gatk.hg38.vcf.gz"),
        variant_caller=variant_caller_gatk)

    bam_file_2 = BamFile(
        path=seq_run_path("1_BAM/fake_sample_2.hg38.bam"),
        aligner=aligner)
    vcf_file_2 = VCFFile(
        path=seq_run_path("2_variants/gatk_per_sample/fake_sample_2.gatk.hg38.vcf.gz"),
        variant_caller=variant_caller_gatk)

    sequencing_files = [
        SequencingFile(sample_name="fake_sample_1",
                       fastq_r1=seq_run_path("0_fastq/fake_sample_1_R1.fastq.gz"),
                       fastq_r2=seq_run_path("0_fastq/fake_sample_1_R2.fastq.gz"),
                       bam_file=bam_file_1,
                       vcf_file=vcf_file_1),
        SequencingFile(sample_name="fake_sample_2",
                       fastq_r1=seq_run_path("0_fastq/fake_sample_2_R1.fastq.gz"),
                       fastq_r2=seq_run_path("0_fastq/fake_sample_2_R1.fastq.gz"),
                       bam_file=bam_file_2,
                       vcf_file=vcf_file_2)
    ]

    qc_by_sample_name = _get_qc_by_sample_name(sample_sheet_lookup, sequencing_files)
    gene_list = [
        "TUBA1A",
        "TUBA8",
        "FLNA",
        "TUBB2B",
        "TUBB3",
        "COL4A1",
        "KIAA1279"
    ]
    qc_gene_lists = [
        QCGeneList(
            path=seq_run_path(f"0_goi/{SEQUENCING_RUN_NAME}_fake_sample_1.txt"),
            qc=qc_by_sample_name["fake_sample_1"],
            gene_list=gene_list),
        QCGeneList(
            path=seq_run_path(f"0_goi/{SEQUENCING_RUN_NAME}_fake_sample_2.txt"),
            qc=qc_by_sample_name["fake_sample_2"],
            gene_list=gene_list)
    ]

    qc_exec_stats_filename_1 = seq_run_path('4_QC/exec_stats/fake_sample_1_qc_summary.txt')
    qc_exec_stats_filename_2 = seq_run_path('4_QC/exec_stats/fake_sample_2_qc_summary.txt')

    qc_exec_stats = [
        QCExecStats(
            qc=qc_by_sample_name["fake_sample_1"],
            path=qc_exec_stats_filename_1, hash=_file_md5sum(qc_exec_stats_filename_1),
            created=datetime.fromisoformat('2024-10-18T11:45:26.826823+10:30'),
            modified=datetime.fromisoformat('2024-10-18T11:45:26.826844+10:30'),
            is_valid=True, deduplicated_reads=165107, indels_dbsnp_percent=95.75, mean_coverage_across_genes=162.84,
            mean_coverage_across_kit=201.43, median_insert=153.0, number_indels=923, number_snps=363,
            percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=98.35, percent_error_rate=0.91,
            percent_map_to_diff_chr=0.75, percent_read_enrichment=54.28, percent_reads=3.552, percent_softclip=0.02,
            percent_duplication=3.42, reads=120760, sample_id_lod=16.6, sex_match='M=yes', snp_dbsnp_percent=96.62,
            ts_to_tv_ratio=2.1, uniformity_of_coverage=84.69),
        QCExecStats(
            qc=qc_by_sample_name["fake_sample_2"],
            path=qc_exec_stats_filename_2, hash=_file_md5sum(qc_exec_stats_filename_2),
            created=datetime.fromisoformat('2024-10-18T11:45:26.838244+10:30'),
            modified=datetime.fromisoformat('2024-10-18T11:45:26.838262+10:30'),
            is_valid=True, deduplicated_reads=275107, indels_dbsnp_percent=88.75, mean_coverage_across_genes=162.84,
            mean_coverage_across_kit=150.43, median_insert=222.0, number_indels=853, number_snps=1213,
            percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=87.35, percent_error_rate=0.55,
            percent_map_to_diff_chr=0.75, percent_read_enrichment=55.28, percent_reads=3.32, percent_softclip=0.02,
            percent_duplication=5.42, reads=410760, sample_id_lod=16.6, sex_match='M=yes', snp_dbsnp_percent=88.62,
            ts_to_tv_ratio=2.1, uniformity_of_coverage=83.69)
    ]

    gene_coverage_filename = seq_run_path("4_QC/bam_stats/samples/fake_sample_1.per_gene_coverage.tsv.gz")

    qc_gene_coverage_list = [
        QCGeneCoverage(qc=qc_by_sample_name["fake_sample_1"],
                       path=gene_coverage_filename),
        QCGeneCoverage(qc=qc_by_sample_name["fake_sample_2"],
                       path=seq_run_path("4_QC/bam_stats/samples/fake_sample_2.per_gene_coverage.tsv.gz"))
    ]

    #########################
    # Call API

    vg_api = VariantGridAPI(server, api_token)

    API_STEPS = {
        "experiment": lambda: vg_api.create_experiment(experiment),
        "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
        "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
        "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
        "sample_sheet_combined_vcf_file": lambda: vg_api.create_sample_sheet_combined_vcf_file(
            sample_sheet_combined_vcf_file),
        "sequencing_data": lambda: vg_api.create_sequencing_data(sample_sheet_lookup, sequencing_files),
        "qc_gene_list": lambda: vg_api.create_qc_gene_list(qc_gene_lists[0]),
        "qc_gene_lists": lambda: vg_api.create_multiple_qc_gene_lists(qc_gene_lists),
        "qc_exec_summary": lambda: vg_api.create_qc_exec_stats(qc_exec_stats[0]),
        "qc_exec_summaries": lambda: vg_api.create_multiple_qc_exec_stats(qc_exec_stats),
        "qc_gene_coverage": lambda: vg_api.create_multiple_qc_gene_coverage(qc_gene_coverage_list),
        "upload_qc_gene_coverage_file": lambda: vg_api.upload_file(gene_coverage_filename),
        "sequencing_run_has_any_vcf": lambda: vg_api.sequencing_run_has_vcf(sequencing_run),
        "sequencing_run_has_fake_vcf": lambda: vg_api.sequencing_run_has_vcf(sequencing_run, "fake_vcf"),
        "sequencing_run_has_our_vcf": lambda: vg_api.sequencing_run_has_vcf(sequencing_run, combo_vcf_filename),
        "upload_vcf_file": lambda: vg_api.upload_file(combo_vcf_filename),
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

