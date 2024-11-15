#!/bin/env python3

import argparse
import os
import re

from samshee.sectionedsheet import read_sectionedsheet
from run_parameters import get_run_parameters
from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import SequencingRun, EnrichmentKit


# Utils

def get_enrichment_kit_and_version_from_scheme(scheme: str) -> tuple[str, int]:
    # scheme looks like 'idt_haem_v1'
    if m := re.match(r"^(.*?)_v(\d)$", scheme):
        ek_name = m.group(1)
        ek_version = int(m.group(2))
    else:
        raise ValueError("Could not extract EnrichmentKit/version from scheme (expected <KIT_NAME>_v<INTEGER>)")
    return ek_name, ek_version


###########
def _get_api(args) -> VariantGridAPI:
    return VariantGridAPI(args.server, args.api_token)

def experiment(args):
    print(f"Running experiment with args: {args}")
    vg_api = _get_api(args)
    vg_api.create_experiment(args.name)

def enrichment_kit(args):
    print(f"Processing enrichment kit with args: {args}")
    vg_api = _get_api(args)
    enrichment_kit = EnrichmentKit(name=args.name,
                                   version=args.version)
    vg_api.create_enrichment_kit(enrichment_kit)


def _get_sequencing_run(sequencing_run_dir, run_parameters_dir, scheme):
    path = args.sequencing_run_dir
    name = os.path.basename(path)
    date = SequencingRun.get_date_from_name(name)

    instrument_name, experiment_name = get_run_parameters(args.run_parameters_dir)

    ek_name, ek_version = get_enrichment_kit_and_version_from_scheme(args.scheme)
    enrichment_kit = EnrichmentKit(name=ek_name,
                                   version=ek_version)
    sequencing_run = SequencingRun(path=path,
                                   name=name,
                                   date=date,
                                   sequencer=instrument_name,
                                   experiment=experiment_name,
                                   enrichment_kit=enrichment_kit)
    return sequencing_run


def _get_sample_sheet(filename):
    sectionedsheet = read_sectionedsheet(filename)
    sequencing_samples = []
    for i, row in enumerate(sectionedsheet["Data"]):
        ss = SequencingSample(sample_id=row["Sample_ID"],
                              sample_project=row["Sample_Project"],
                              sample_number=i+1,
                              lane=row.get("Lane"),
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


def sequencing_run(args):
    print(f"Running sequencing run with args: {args}")
    # TAU seem to pass in enrichment kit version as follows
    # --scheme idt_haem_v1

    sequencing_run = _get_sequencing_run(args.sequencing_run_dir, args.run_parameters_dir, args.scheme)
    vg_api = _get_api(args)
    vg_api.create_sequencing_run(sequencing_run)

def sample_sheet(args):
    print(f"Processing sample sheet with args: {args}")

def sample_sheet_combined_vcf_file(args):
    print(f"Handling combined VCF file for sample sheet with args: {args}")

def sequencing_data(args):
    print(f"Handling sequencing data with args: {args}")

def qc_gene_list(args):
    print(f"Processing QC gene_list list with args: {args}")

def qc_exec_stats(args):
    print(f"Executing QC stats with args: {args}")

def qc_gene_coverage(args):
    print(f"Checking QC gene coverage with args: {args}")

def upload_file(args):
    print(f"Uploading file with args: {args}")

def sequencing_run_all_files(args):
    pass



def parse_args():
    parser = argparse.ArgumentParser(description="Bioinformatics Command Line Interface")
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--api-token', required=True, help='API token for authentication')


    subparsers = parser.add_subparsers(help="Available subcommands", required=True)
    # EVERYTHING
    subparser_experiment = subparsers.add_parser("sequencing_run_all_files",
                                                 help="")
    subparser_experiment.set_defaults(func=sequencing_run_all_files)
    subparser_experiment.add_argument("--path", help="SequencingRun dir", required=True)
    subparser_experiment.add_argument("--run-parameters-dir", help="Dir containing RunParameters", required=True)
    subparser_experiment.add_argument("--scheme", help="enrichment kit and version", required=True)

    # Experiment
    subparser_experiment = subparsers.add_parser("experiment",
                                                 help="Manually creating experiment (this is usually done as part of 'sequencing_run')")
    subparser_experiment.set_defaults(func=experiment)
    subparser_experiment.add_argument("--name", help="Name of experiment to create", required=True)

    # EnrichmentKit
    subparser_experiment = subparsers.add_parser("enrichment_kit",
                                                 help="Manually creating enrichment_kit (this is usually done as part of 'sequencing_run')")
    subparser_experiment.set_defaults(func=enrichment_kit)
    subparser_experiment.add_argument("--name", help="Name of enrichment_kit", required=True)
    subparser_experiment.add_argument("--version", help="version", default=10, type=int)

    # SequencingRun
    subparser_experiment = subparsers.add_parser("sequencing_run",
                                                 help="")
    subparser_experiment.set_defaults(func=sequencing_run)
    subparser_experiment.add_argument("--path", help="SequencingRun dir", required=True)
    subparser_experiment.add_argument("--run-parameters-dir", help="Dir containing RunParameters", required=True)
    subparser_experiment.add_argument("--scheme", help="enrichment kit and version", required=True)

    # SampleSheet
    subparser_experiment = subparsers.add_parser("sample_sheet",
                                                 help="")
    subparser_experiment.set_defaults(func=sample_sheet)

    # SampleSheet
    subparser_experiment = subparsers.add_parser("sample_sheet_combined_vcf_file",
                                                 help="")
    subparser_experiment.set_defaults(func=sample_sheet_combined_vcf_file)

    # sequencing_data
    subparser_experiment = subparsers.add_parser("sequencing_data",
                                                 help="")
    subparser_experiment.set_defaults(func=sequencing_data)

    # qc_gene_list
    subparser_experiment = subparsers.add_parser("qc_gene_list",
                                                 help="")
    subparser_experiment.set_defaults(func=qc_gene_list)

    # qc_exec_stats
    subparser_experiment = subparsers.add_parser("qc_exec_stats",
                                                 help="")
    subparser_experiment.set_defaults(func=qc_exec_stats)

    # qc_gene_coverage
    subparser_experiment = subparsers.add_parser("qc_gene_coverage",
                                                 help="")
    subparser_experiment.set_defaults(func=qc_gene_coverage)


    # upload_file
    subparser_experiment = subparsers.add_parser("upload_file",
                                                 help="")
    subparser_experiment.set_defaults(func=upload_file)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    args.func(args)
