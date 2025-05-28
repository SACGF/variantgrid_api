import argparse
import dataclasses
import os
from hashlib import md5
from typing import List, Dict
import re
from datetime import datetime
import pandas as pd
import csv


from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit, SequencingRun, SequencingSample, SampleSheet, \
    SampleSheetCombinedVCFFile, VariantCaller, SampleSheetLookup, Aligner, VCFFile, BamFile, SequencingFile, \
    SequencingSampleLookup, QC, QCGeneList, QCExecStats, QCGeneCoverage


def addArgs(parser):
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--sample', required=True, help='Sample ID')
    #parser.add_argument('--sample_sheet', required=True, help='Samplesheet')
    parser.add_argument('--api-token', required=True, help='API token for authentication')
    parser.add_argument('--step', required=False, help='Run a single step (default: run all)')
    parser.add_argument('--batch_id', help= 'batch ID of the run (LabrunID + Sequencing run ID)', required=True)
    parser.add_argument('--results_dir', help= 'Final destination path of run results folder', required=True)
    parser.add_argument('--continue_on_error', help= 'Optionally continue even if errors; default dies immediately', action='store_true')
    parser.add_argument('--api_complete', help= 'Create final output once API is finished', required=False)
    parser.add_argument('--api_log', help= 'If specified, output log file', required=False)

    args = parser.parse_args()
    return args

parser = argparse.ArgumentParser(description="VariantGrid API client")
args = addArgs(parser)

# get arguments and define variables

BATCHID = args.batch_id
API_TOKEN = args.api_token
STEP = args.step
SERVER = args.server
TAU_RESULTS_DIR = args.results_dir
SAMPLE = args.sample
#SAMPLESHEET =args.sample_sheet
API_COMP = args.api_complete
API_LOG = args.api_log
CONTINUE_ON_ERROR = args.continue_on_error


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





def get_sample_info_from_samplesheet(samplesheet_path, sample_id):

    # Open the file and find the line where "FlowCellID" appears
    with open(samplesheet_path, "r") as f:
        lines = f.readlines()

    # Find the index of the line containing "FlowCellID"
    start_index = None
    for i, line in enumerate(lines):
        if "FlowCellID" in line:
            start_index = i
            break

    if start_index is None:
        raise ValueError("FlowCellID not found in the SampleSheet.csv file.")

    # Read the file into a DataFrame starting from the "FloCellID" line
    df = pd.read_csv(samplesheet_path, skiprows=start_index)

    # Search for the SAMPLE in the Sample_ID column
    sample_row = df[df["Sample_ID"] == sample_id]

    if sample_row.empty:
        raise ValueError(f"Sample_ID '{sample_id}' not found in the SampleSheet.csv file.")

    # Retrieve the row number and values of two other columns
    row_number = int(sample_row.index[0])
    sapnumber = sample_row.iloc[0]["SAPOrderNumber"]  
    sample_index = sample_row.iloc[0]["index"]

    # Set control to True or False based on sapnumber
    control = sapnumber == "Control"

    return row_number,  sapnumber, sample_index,control





def get_experiment_info(BATCHID):
    base_name = BATCHID

    parts = base_name.split("_")

    experiment = "_".join(parts[:3])   

    date_obj = datetime.strptime(parts[3], "%y%m%d")
    formatted_date = date_obj.strftime("%Y-%m-%d")  

    sequencer = parts[4]

    #kit = os.path.basename(os.path.dirname(TAU_RESULTS_DIR))

    ##Kit and version info, if new kit added, needs to be added here. TWIST coming soon....
    if('Exome' in experiment):
            kit_name = 'idt_exome'
            kit_name_version = 2

    elif('GMPFocus' in experiment):
        kit_name = 'idt_gmp_focus'
        kit_name_version = 4_1

    elif('BRCAcan' in experiment):
        kit_name = 'idt_brca_cancer'
        kit_name_version = 1


    elif('Haem' in experiment):
        kit_name = 'idt_haem'
        kit_name_version = 1
    elif('RhampFFPEonco' in experiment):
        kit_name = 'idt_rhampseq_ffpe_onco'
        kit_name_version = 2
    elif('RhampIDHplus' in experiment):
        kit_name = 'idt_rhampseq_idhplus'
        kit_name_version = 2


    return experiment.upper(), formatted_date, sequencer, kit_name, kit_name_version


def parse_exec_stats(tsv_path):

    stats = {}
    with open(tsv_path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) >= 2:
                stats[row[0].strip()] = row[1].strip()
    return stats




def run_api(server, api_token, step=None):
    


    #samplesheet_file = f'{seq2_run_path}/SampleSheet.csv'

    #data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
    seq_run_dir = os.path.join(TAU_RESULTS_DIR, BATCHID)
    def seq_run_path(path):
       return os.path.join(seq_run_dir, path)
    #Get run info from the batch ID
    experiment, formatted_date, sequencer, kit_name, kit_name_version = get_experiment_info(BATCHID)

    samplesheet_path = seq_run_path("SampleSheet.csv")
    sample_id = SAMPLE  # Replace SAMPLE with the actual sample ID

    # Get sample info from the samplesheet
    row_number,  sapnumber, sample_index,control = get_sample_info_from_samplesheet(samplesheet_path, sample_id)


    date_obj = datetime.strptime(formatted_date, "%Y-%m-%d")
    unix_timestamp = date_obj.timestamp()






    #experiment = "HAEM_20_999"
    enrichment_kit = EnrichmentKit(name=kit_name, version=kit_name_version)
    sequencing_run = SequencingRun(path=seq_run_dir,
                                   name=BATCHID,
                                   date=formatted_date,
                                   sequencer=sequencer,
                                   experiment=experiment,
                                   enrichment_kit=enrichment_kit)

    sequencing_samples = [
        SequencingSample(sample_id=SAMPLE,
                         sample_project=None,
                         sample_number=row_number,
                         #lane=3, #optional field
                         barcode=sample_index, 
                         enrichment_kit=enrichment_kit,
                         is_control=control,
                         failed=False,
                         data=[
                             {
                                 "column": "SAPOrderNumber",
                                 "value": sapnumber
                             }
                         ])
    ]

    sample_sheet = SampleSheet(
        path=seq_run_path("SampleSheet.csv"),
        sequencing_run=sequencing_run,
        file_last_modified=unix_timestamp,
        hash=_file_md5sum(seq_run_path("SampleSheet.csv")),#"f0ac87bcae3f0e56b3f65b70fd6389ce"
        sequencing_samples=sequencing_samples)

    sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)

    variant_caller_var_dict = VariantCaller(name="VarDict", version="1.8.2")
    combo_vcf_filename = seq_run_path(f"2_variants/{BATCHID}.vardict.hg38.vcf.gz")
    sample_sheet_combined_vcf_file = SampleSheetCombinedVCFFile(
        path=combo_vcf_filename,
        sample_sheet_lookup=sample_sheet_lookup,
        variant_caller=variant_caller_var_dict)

    aligner = Aligner(name='BWA', version="0.7.18")
    variant_caller_gatk = VariantCaller(name="GATK", version="4.1.9.0")
    bam_file = BamFile(
        path=seq_run_path(f"1_BAM/{SAMPLE}.hg38.bam"),
        aligner=aligner)
    vcf_file = VCFFile(
        path=seq_run_path(f"2_variants/gatk_per_sample/{SAMPLE}.gatk.hg38.vcf.gz"),
        variant_caller=variant_caller_gatk)

    sequencing_files = [
        SequencingFile(sample_name=f"{SAMPLE}",
                       fastq_r1=seq_run_path(f"0_fastq/{SAMPLE}_R1.fastq.gz"),
                       fastq_r2=seq_run_path(f"0_fastq/{SAMPLE}_R2.fastq.gz"),
                       bam_file=bam_file,
                       vcf_file=vcf_file)

    ]

    qc_by_sample_name = _get_qc_by_sample_name(sample_sheet_lookup, sequencing_files)

    
    


    ##This will use goi in the goi folder and if not present, will use the one defined in the enrichment kit
    gene_list_path = seq_run_path(f"0_goi/{BATCHID}_{SAMPLE}.txt")
    
    if os.path.exists(gene_list_path):
        with open(gene_list_path) as f:
            gene_list = [line.strip() for line in f if line.strip()]
    else:
        vg_api = VariantGridAPI(server, api_token)
        enrichment_kit_response = vg_api.create_enrichment_kit(enrichment_kit)
        gene_list = [
            g['gene_symbol']
            for g in enrichment_kit_response['gene_list']['genelistgenesymbol_set']
            if 'gene_symbol' in g and g['gene_symbol']
        ]
    qc_gene_lists = [
        QCGeneList(
            path=seq_run_path(f"0_goi/{BATCHID}_{SAMPLE}.txt"),
            qc=qc_by_sample_name[f"{SAMPLE}"],
            gene_list=gene_list)
    ]

    qc_exec_stats_file = seq_run_path(f"4_QC/exec_stats/{SAMPLE}_qc_summary.txt")

    sample_exec_stats = seq_run_path(f"4_QC/exec_stats/{SAMPLE}_stats.tsv")

    exec_stats_dict = parse_exec_stats(sample_exec_stats)


    if BATCHID.startswith("Exome") or BATCHID.startswith("GMPFocus"):
        qc_exec_stats = [
            QCExecStats(
                qc=qc_by_sample_name[f"{SAMPLE}"],
                path=qc_exec_stats_file, hash=_file_md5sum(qc_exec_stats_file),
                created='2024-10-18T11:45:26.826823+10:30', modified='2024-10-18T11:45:26.826844+10:30',
                is_valid=True,
                deduplicated_reads=165107,
                indels_dbsnp_percent=95.75, 
                mean_coverage_across_genes=162.84,
                mean_coverage_across_kit=201.43,
                median_insert=153.0,
                number_indels=923,
                number_snps=363,
                percent_10x_goi=100.0,
                percent_20x_goi=100.0,
                percent_20x_kit=98.35,
                percent_error_rate=0.91,
                percent_map_to_diff_chr=0.75,
                percent_read_enrichment=54.28,
                percent_reads=3.552,percent_softclip=0.02,
                percent_duplication=3.42,
                reads=120760,
                sample_id_lod=16.6,
                sex_match='M=yes', 
                snp_dbsnp_percent=96.62,
                ts_to_tv_ratio=2.1,
                uniformity_of_coverage=84.69
            )
        ]   
    
    elif BATCHID.startswith("RhampFFPE") or BATCHID.startswith("RhampIDHplus") or BATCHID.startswith("Haem"):
        qc_exec_stats = [
            QCExecStats(
                qc=qc_by_sample_name[f"{SAMPLE}"],
                path=qc_exec_stats_file, hash=_file_md5sum(qc_exec_stats_file),
                created='2024-10-18T11:45:26.826823+10:30', modified='2024-10-18T11:45:26.826844+10:30',
                is_valid=True,
                deduplicated_reads=int(exec_stats_dict.get('DedupReads', 0)),
                indels_dbsnp_percent=0.0,
                mean_coverage_across_genes=float(exec_stats_dict.get('MeanDepth_GOI', 0)),
                mean_coverage_across_kit=float(exec_stats_dict.get('MeanDepth_Kit', 0)),
                median_insert=float(exec_stats_dict.get('MedianInsert', 0)),
                number_indels=0,
                number_snps=0,
                percent_10x_goi=100.0,
                percent_20x_goi=100.0,
                percent_20x_kit=98.35,
                percent_error_rate=0.91,
                percent_map_to_diff_chr=0.75,
                percent_read_enrichment=54.28,
                percent_reads=3.552,
                percent_softclip=0.02,
                percent_duplication=3.42,
                reads=120760,
                sample_id_lod=16.6,
            )
        ]
    gene_coverage_filename = seq_run_path(f"4_QC/bam_stats/samples/{SAMPLE}.per_gene_coverage.tsv.gz")
    print(f"gene_coverage_filename: {gene_coverage_filename}")
    qc_gene_coverage_list = [
        QCGeneCoverage(qc=qc_by_sample_name[f"{SAMPLE}"],
                       path=gene_coverage_filename)
    ]

    #########################
    # Call API

    vg_api = VariantGridAPI(server, api_token)

    import dataclasses
    #print("dict here ", dataclasses.asdict(sample_sheet))

    #print(f"{combo_vcf_filename=}")

    API_STEPS = {
        "experiment": lambda: vg_api.create_experiment(experiment),
        "enrichment_kit": lambda: enrichment_kit_response,
        "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
        "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
        "sample_sheet_combined_vcf_file": lambda: vg_api.create_sample_sheet_combined_vcf_file(
            sample_sheet_combined_vcf_file),
        "sequencing_data": lambda: vg_api.create_sequencing_data(sample_sheet_lookup, sequencing_files),
        "qc_gene_list": lambda: vg_api.create_qc_gene_list(qc_gene_lists[0]),
        #"qc_gene_lists": lambda: vg_api.create_multiple_qc_gene_lists(qc_gene_lists), ##creates multiple gene lists, we don't need this since we are doing one sample at a time
        "qc_exec_summary": lambda: vg_api.create_qc_exec_stats(qc_exec_stats[0]),
        "qc_exec_summaries": lambda: vg_api.create_multiple_qc_exec_stats(qc_exec_stats),
        "qc_gene_coverage": lambda: vg_api.create_multiple_qc_gene_coverage(qc_gene_coverage_list),
        "upload_qc_gene_coverage_file": lambda: vg_api.upload_file(gene_coverage_filename),
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
    run_api(SERVER, API_TOKEN,
             step=args.step)

