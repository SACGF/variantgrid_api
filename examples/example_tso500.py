"""
    A TSO 500 run: one specimen sequenced as a DNA arm and an RNA arm.

    1. Accession the patient, the specimen, and its two extractions (DNA '2600000001C', RNA '2600000001B')
    2. Post the sequencing run and sample sheet as usual
    3. Link each arm's sequencing sample to its extraction - one call per arm
    4. Upload the DNA arm's VCFs and the RNA arm's files, naming the extraction and build in upload metadata
    5. Post the specimen's TMB / MSI / GIS, transcribed from CombinedVariantOutput.tsv

    Ordering is forgiving: a link or upload naming an extraction the server doesn't have yet is parked and
    attaches itself once the extraction is created, so nothing needs re-sending.

    Needs a server at or after SACGF/variantgrid#1716. Test data is synthetic - see tests/test_data/tso500/README.md
"""
import argparse
import os
from datetime import datetime
from typing import Dict

from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit, SequencerModel, Sequencer, SequencingRun, SequencingSample, \
    SampleSheet, SampleSheetLookup, SequencingSampleLookup, Patient, Specimen, Extraction, \
    SpecimenMeasure, ExternalReference, TissueStatus, NucleicAcid, SpecimenMeasureType


def parse_args():
    parser = argparse.ArgumentParser(description="VariantGrid API client - TSO 500")
    parser.add_argument('--server', required=True, help='Base URL of the VariantGrid server inc. port')
    parser.add_argument('--api-token', required=True, help='API token for authentication')
    parser.add_argument('--step', required=False, help='Run a single step (default: run all)')
    return parser.parse_args()


def read_combined_variant_output(filename) -> Dict[str, Dict[str, str]]:
    """ {section: {key: value}} for the key/value sections of a CombinedVariantOutput.tsv, eg [TMB] """
    sections = {}
    section = None
    with open(filename) as f:
        for line in f:
            columns = [c for c in line.rstrip("\n").split("\t") if c]
            if not columns:
                continue
            if columns[0].startswith("[") and columns[0].endswith("]"):
                section = sections.setdefault(columns[0][1:-1], {})
            elif section is not None and len(columns) == 2:
                section[columns[0]] = columns[1]
    return sections


def get_specimen_measures(cvo_filename, dna_extraction: ExternalReference):
    """ Transcribe, never compute: the values are copied out of the vendor file with the block they came
        from as source_payload. Set call / threshold / threshold_source only from your lab's own policy. """
    sections = read_combined_variant_output(cvo_filename)
    analysis_details = sections["Analysis Details"]
    method = f"DRAGEN TSO500 CombinedVariantOutput {analysis_details['Module Version']}"
    measured_date = datetime.fromisoformat(f"{analysis_details['Output Date']}T{analysis_details['Output Time']}")

    # (measure_type, section, key, unit)
    measure_sources = [
        (SpecimenMeasureType.TMB, "TMB", "Total TMB", "mut/Mb"),
        (SpecimenMeasureType.MSI, "MSI", "Percent Unstable MSI Sites", "%"),
        (SpecimenMeasureType.GIS, "GIS", "Genomic Instability Score", None),
        (SpecimenMeasureType.TUMOUR_FRACTION, "GIS", "Tumor Fraction", None),
        (SpecimenMeasureType.PLOIDY, "GIS", "Ploidy", None),
    ]
    measures = []
    for measure_type, section, key, unit in measure_sources:
        values = sections.get(section, {})
        if values.get(key) in (None, "NA"):
            continue
        measures.append(SpecimenMeasure(measure_type=measure_type,
                                        value=float(values[key]),
                                        unit=unit,
                                        method=method,
                                        source_payload=values,
                                        measured_date=measured_date,
                                        extraction=dna_extraction))
    return measures


def test_api(server, api_token, step=None):
    proj_dir = os.path.dirname(os.path.dirname(__file__))  # In "examples"
    data_dir = os.path.join(proj_dir, "tests", "test_data", "tso500", "ExampleSample_2600000001")
    dna_dir = os.path.join(data_dir, "ExampleSample_DNA_2600000001C")
    rna_dir = os.path.join(data_dir, "ExampleSample_RNA_2600000001B")

    # The pair ID is the patient, the 10 digit accession the specimen, and the C / B container suffixes
    # the specimen's DNA and RNA extractions
    patient = Patient(patient_code="C0000001", affected=True)
    specimen = Specimen(patient=ExternalReference.from_patient(patient),
                        reference_id="2600000001",
                        tissue_status=TissueStatus.AFFECTED)
    specimen_reference = ExternalReference.from_specimen(specimen)
    dna_extraction = Extraction(specimen=specimen_reference, reference_id="2600000001C",
                                nucleic_acid_source=NucleicAcid.DNA)
    rna_extraction = Extraction(specimen=specimen_reference, reference_id="2600000001B",
                                nucleic_acid_source=NucleicAcid.RNA)
    dna_reference = ExternalReference.from_extraction(dna_extraction)
    rna_reference = ExternalReference.from_extraction(rna_extraction)

    # Sequencing run - both arms are on the one sample sheet
    run_name = "TSO500_26_001_260101_LH00639_0100_A23ABCDEF1"
    enrichment_kit = EnrichmentKit(name="tso500", version=1)
    sequencer = Sequencer(name="LH00639", sequencer_model=SequencerModel.from_sequencer_name("LH00639"))
    sequencing_run = SequencingRun(path=os.path.join("/data/tso500", run_name),
                                   name=run_name,
                                   date=SequencingRun.get_date_from_name(run_name),
                                   sequencer=sequencer.name,
                                   experiment="TSO500_26_001",
                                   enrichment_kit=enrichment_kit)
    dna_sample_name = "ExampleSample_DNA_2600000001C"
    rna_sample_name = "ExampleSample_RNA_2600000001B"
    sample_sheet = SampleSheet(path=os.path.join(sequencing_run.path, "SampleSheet.csv"),
                               sequencing_run=sequencing_run,
                               file_last_modified=1767225600.0,
                               hash="5d1e4c6f0b2a8e7d9c3b1a0f6e5d4c3b",
                               sequencing_samples=[
                                   SequencingSample(sample_id=dna_sample_name, sample_number=1, barcode="GCCAAT",
                                                    enrichment_kit=enrichment_kit),
                                   SequencingSample(sample_id=rna_sample_name, sample_number=2, barcode="CAGATC",
                                                    enrichment_kit=enrichment_kit),
                               ])
    sample_sheet_lookup = SampleSheetLookup.from_sample_sheet(sample_sheet)
    dna_sequencing_sample = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup,
                                                   sample_name=dna_sample_name)
    rna_sequencing_sample = SequencingSampleLookup(sample_sheet_lookup=sample_sheet_lookup,
                                                   sample_name=rna_sample_name)

    # Upload metadata - facts the files don't reliably carry. Send a build's own name (GRCh37), not an alias (hg19)
    dna_prefix = os.path.join(dna_dir, dna_sample_name)
    rna_prefix = os.path.join(rna_dir, rna_sample_name)
    uploads = {
        "small_variants": (f"{dna_prefix}.hard-filtered.vcf",
                           {"extraction": dna_reference, "source": "DRAGEN TSO500 SmallVariant"}),
        "cnv": (f"{dna_prefix}.cnv.vcf",
                {"extraction": dna_reference, "source": "DRAGEN TSO500 CNV"}),
        # No contigs in the header, so the build must be declared
        "exon_cnv": (f"{dna_prefix}_DragenExonCNV.vcf",
                     {"extraction": dna_reference, "genome_build": "GRCh37"}),
        "splice_variants": (f"{rna_prefix}_SpliceVariants.vcf",
                            {"extraction": rna_reference}),
        # Carries no build at all
        "fusions": (f"{rna_prefix}_AllFusions.csv",
                    {"extraction": rna_reference, "genome_build": "GRCh37"}),
    }

    cvo_filename = os.path.join(data_dir, "ExampleSample_2600000001_CombinedVariantOutput.tsv")
    specimen_measures = get_specimen_measures(cvo_filename, dna_reference)

    #########################
    # Call API

    vg_api = VariantGridAPI(server, api_token)

    API_STEPS = {
        # 1. Accessioning - a specimen needs its patient, an extraction its specimen
        "patient": lambda: vg_api.create_patient(patient),
        "specimen": lambda: vg_api.create_specimen(specimen),
        "dna_extraction": lambda: vg_api.create_extraction(dna_extraction),
        "rna_extraction": lambda: vg_api.create_extraction(rna_extraction),
        # 2. Sequencing
        "enrichment_kit": lambda: vg_api.create_enrichment_kit(enrichment_kit),
        "sequencer_model": lambda: vg_api.create_sequencer_model(sequencer.sequencer_model),
        "sequencer": lambda: vg_api.create_sequencer(sequencer),
        "sequencing_run": lambda: vg_api.create_sequencing_run(sequencing_run),
        "sample_sheet": lambda: vg_api.create_sample_sheet(sample_sheet),
        # 3. One link per arm. Returns match_status 'Pending' (HTTP 202) if the extraction isn't there yet
        "link_dna_extraction": lambda: vg_api.link_sequencing_sample_extraction(dna_sequencing_sample, dna_reference),
        "link_rna_extraction": lambda: vg_api.link_sequencing_sample_extraction(rna_sequencing_sample, rna_reference),
    }
    # 4. Uploads. path=None as these aren't registered SeqAuto VCFs - the metadata names the extraction instead
    for name, (filename, metadata) in uploads.items():
        API_STEPS[f"upload_{name}"] = lambda f=filename, m=metadata: vg_api.upload_file(f, path=None, metadata=m)
    # 5. Measures describe the specimen; the DNA arm produced them
    API_STEPS["specimen_measures"] = lambda: vg_api.create_specimen_measures(specimen_reference, specimen_measures)

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
