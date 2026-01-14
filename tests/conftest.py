import os
from datetime import datetime
from pathlib import Path
import pytest

from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import (
    EnrichmentKit, SequencerModel, Sequencer, SequencingRun, SequencingSample, SampleSheet,
    SampleSheetCombinedVCFFile, VariantCaller, SampleSheetLookup, Aligner, VCFFile, BamFile,
    SequencingFile, SequencingSampleLookup, QC, QCGeneList, QCExecStats, QCGeneCoverage, Manufacturer
)

@pytest.fixture
def server(): return "https://example.org"

@pytest.fixture
def api_token(): return "TKN"

@pytest.fixture
def api(server, api_token): return VariantGridAPI(server, api_token)

@pytest.fixture
def data_dir(): return Path(__file__).parent / "test_data"

@pytest.fixture
def vg_objects(data_dir):
    SEQUENCING_RUN_NAME = "Haem_20_999_201231_M02027_0112_000000000_JFT79"
    seq_run_dir = data_dir / "idt_haem" / SEQUENCING_RUN_NAME
    seq_run_path = lambda p: str(seq_run_dir / p)

    experiment = "HAEM_20_999"
    enrichment_kit = EnrichmentKit(name="idt_haem", version=1)

    manufacturer = Manufacturer(name="Illumina")
    sequencer_model = SequencerModel(model="HiSeq 2500", manufacturer=manufacturer, data_naming_convention="H")
    sequencer = Sequencer(name="SN1101", sequencer_model=sequencer_model)

    seq_date = SequencingRun.get_date_from_name(SEQUENCING_RUN_NAME)
    sequencing_run = SequencingRun(
        path=str(seq_run_dir), name=SEQUENCING_RUN_NAME, date=seq_date, sequencer=sequencer.name,
        experiment=experiment, enrichment_kit=enrichment_kit
    )

    sequencing_samples = [
        SequencingSample(sample_id="fake_sample_1", sample_project=None, sample_number=1, barcode="GCCAAT",
                         enrichment_kit=enrichment_kit, is_control=False, failed=False,
                         data=[{"column": "SAPOrderNumber", "value": "SAP1000001"}]),
        SequencingSample(sample_id="fake_sample_2", sample_project=None, sample_number=2, barcode="CAGATC",
                         enrichment_kit=enrichment_kit, is_control=False, failed=False,
                         data=[{"column": "SAPOrderNumber", "value": "SAP1000002"}]),
    ]

    sample_sheet = SampleSheet(
        path=seq_run_path("SampleSheet.csv"),
        sequencing_run=sequencing_run,
        file_last_modified=1725941707.0033002,
        hash="f0ac87bcae3f0e56b3f65b70fd6389ce",
        sequencing_samples=sequencing_samples
    )

    sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)

    variant_caller_var_dict = VariantCaller(name="VarDict", version="1.8.2")
    combo_vcf_filename = seq_run_path(f"2_variants/{SEQUENCING_RUN_NAME}.vardict.hg38.vcf.gz")
    sample_sheet_combined_vcf_file = SampleSheetCombinedVCFFile(
        path=combo_vcf_filename, sample_sheet_lookup=sample_sheet_lookup, variant_caller=variant_caller_var_dict
    )

    aligner = Aligner(name="BWA", version="0.7.18")
    variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")

    single_sample_vcf_filename_1 = seq_run_path("2_variants/gatk_per_sample/fake_sample_1.gatk.hg38.vcf.gz")
    single_sample_vcf_filename_2 = seq_run_path("2_variants/gatk_per_sample/fake_sample_2.gatk.hg38.vcf.gz")

    bam_file_1 = BamFile(path=seq_run_path("1_BAM/fake_sample_1.hg38.bam"), aligner=aligner)
    vcf_file_1 = VCFFile(path=single_sample_vcf_filename_1, variant_caller=variant_caller_gatk)

    bam_file_2 = BamFile(path=seq_run_path("1_BAM/fake_sample_2.hg38.bam"), aligner=aligner)
    vcf_file_2 = VCFFile(path=single_sample_vcf_filename_2, variant_caller=variant_caller_gatk)

    sequencing_files = [
        SequencingFile(sample_name="fake_sample_1",
                       fastq_r1=seq_run_path("0_fastq/fake_sample_1_R1.fastq.gz"),
                       fastq_r2=seq_run_path("0_fastq/fake_sample_1_R2.fastq.gz"),
                       bam_file=bam_file_1, vcf_file=vcf_file_1),
        SequencingFile(sample_name="fake_sample_2",
                       fastq_r1=seq_run_path("0_fastq/fake_sample_2_R1.fastq.gz"),
                       fastq_r2=seq_run_path("0_fastq/fake_sample_2_R1.fastq.gz"),
                       bam_file=bam_file_2, vcf_file=vcf_file_2),
    ]

    bam_and_vcf = {}
    for sf in sequencing_files:
        bam_and_vcf[sf.sample_name] = (sf.bam_file.__class__(**{**sf.bam_file.__dict__, "aligner": None}),
                                       sf.vcf_file.__class__(**{**sf.vcf_file.__dict__, "variant_caller": None}))

    qc_by_name = {}
    for sample_name, (bam_file, vcf_file) in bam_and_vcf.items():
        ssl = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup, sample_name=sample_name)
        qc_by_name[sample_name] = QC(sequencing_sample_lookup=ssl, bam_file=bam_file, vcf_file=vcf_file)

    gene_list = ["TUBA1A","TUBA8","FLNA","TUBB2B","TUBB3","COL4A1","KIAA1279"]
    qc_gene_lists = [
        QCGeneList(path=seq_run_path(f"0_goi/{SEQUENCING_RUN_NAME}_fake_sample_1.txt"),
                   qc=qc_by_name["fake_sample_1"], gene_list=gene_list),
        QCGeneList(path=seq_run_path(f"0_goi/{SEQUENCING_RUN_NAME}_fake_sample_2.txt"),
                   qc=qc_by_name["fake_sample_2"], gene_list=gene_list),
    ]

    qc_exec_stats = [
        QCExecStats(
            qc=qc_by_name["fake_sample_1"],
            path=seq_run_path("4_QC/exec_stats/fake_sample_1_qc_summary.txt"),
            hash="x",
            created=datetime.fromisoformat("2024-10-18T11:45:26.826823+10:30"),
            modified=datetime.fromisoformat("2024-10-18T11:45:26.826844+10:30"),
            is_valid=True, deduplicated_reads=165107, indels_dbsnp_percent=95.75, mean_coverage_across_genes=162.84,
            mean_coverage_across_kit=201.43, median_insert=153.0, number_indels=923, number_snps=363,
            percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=98.35, percent_error_rate=0.91,
            percent_map_to_diff_chr=0.75, percent_read_enrichment=54.28, percent_reads=3.552, percent_softclip=0.02,
            percent_duplication=3.42, reads=120760, sample_id_lod=16.6, sex_match="M=yes", snp_dbsnp_percent=96.62,
            ts_to_tv_ratio=2.1, uniformity_of_coverage=84.69
        ),
        QCExecStats(
            qc=qc_by_name["fake_sample_2"],
            path=seq_run_path("4_QC/exec_stats/fake_sample_2_qc_summary.txt"),
            hash="y",
            created=datetime.fromisoformat("2024-10-18T11:45:26.838244+10:30"),
            modified=datetime.fromisoformat("2024-10-18T11:45:26.838262+10:30"),
            is_valid=True, deduplicated_reads=275107, indels_dbsnp_percent=88.75, mean_coverage_across_genes=162.84,
            mean_coverage_across_kit=150.43, median_insert=222.0, number_indels=853, number_snps=1213,
            percent_10x_goi=100.0, percent_20x_goi=100.0, percent_20x_kit=87.35, percent_error_rate=0.55,
            percent_map_to_diff_chr=0.75, percent_read_enrichment=55.28, percent_reads=3.32, percent_softclip=0.02,
            percent_duplication=5.42, reads=410760, sample_id_lod=16.6, sex_match="M=yes", snp_dbsnp_percent=88.62,
            ts_to_tv_ratio=2.1, uniformity_of_coverage=83.69
        )
    ]

    qc_gene_coverage_list = [
        QCGeneCoverage(qc=qc_by_name["fake_sample_1"],
                       path=seq_run_path("4_QC/bam_stats/samples/fake_sample_1.per_gene_coverage.tsv.gz")),
        QCGeneCoverage(qc=qc_by_name["fake_sample_2"],
                       path=seq_run_path("4_QC/bam_stats/samples/fake_sample_2.per_gene_coverage.tsv.gz")),
    ]

    return dict(
        SEQUENCING_RUN_NAME=SEQUENCING_RUN_NAME,
        seq_run_dir=str(seq_run_dir),
        experiment=experiment,
        enrichment_kit=enrichment_kit,
        sequencer_model=sequencer_model,
        sequencer=sequencer,
        sequencing_run=sequencing_run,
        sample_sheet=sample_sheet,
        sample_sheet_lookup=sample_sheet_lookup,
        sample_sheet_combined_vcf_file=sample_sheet_combined_vcf_file,
        sequencing_files=sequencing_files,
        qc_gene_lists=qc_gene_lists,
        qc_exec_stats=qc_exec_stats,
        qc_gene_coverage_list=qc_gene_coverage_list,
        combo_vcf_filename=combo_vcf_filename,
        single_sample_vcf_filename_1=single_sample_vcf_filename_1,
        single_sample_vcf_filename_2=single_sample_vcf_filename_2,
    )
