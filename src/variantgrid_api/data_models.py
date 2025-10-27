import re
from datetime import date, datetime
from dataclasses import dataclass, field
from typing import Optional, List

from dataclasses_json import dataclass_json, config


@dataclass_json
@dataclass
class EnrichmentKit:
    name: str
    version: int


@dataclass_json
@dataclass
class SequencerModel:
    model: str
    manufacturer: str
    data_naming_convention: str # 'M' or 'H' for (HiSeq only)


@dataclass_json
@dataclass
class Sequencer:
    name: str
    sequencer_model: SequencerModel


@dataclass_json
@dataclass
class SequencingRun:
    path: str
    name: str
    date: date
    sequencer: str
    experiment: str
    enrichment_kit: EnrichmentKit

    @staticmethod
    def get_date_from_name(name) -> Optional[datetime.date]:
        date = None
        if m := re.match(r".*_?([12]\d{5})_", name):
            date_str = m.group(1)
            dt = datetime.strptime(date_str, "%y%m%d")
            date = dt.date()
        return date



@dataclass_json
@dataclass
class SequencingSample:
    sample_id: str
    sample_number: int  # Sample sheet row (to keep in order)
    barcode: str
    enrichment_kit: EnrichmentKit
    lane: Optional[int] = None
    sample_project: Optional[str] = None
    is_control: bool = False
    failed: bool = False
    data: List[dict] = field(default_factory=lambda: [], metadata=config(field_name="sequencingsampledata_set"))


@dataclass_json
@dataclass
class SampleSheet:
    path: str
    sequencing_run: SequencingRun
    file_last_modified: float
    hash: str
    sequencing_samples: List[SequencingSample] = field(metadata=config(field_name="sequencingsample_set"))


@dataclass_json
@dataclass
class SampleSheetLookup:
    """ 'Lookups' are used as arguments to find existing data on server - not enough details to create one,
        but just enough to find one
    """
    sequencing_run: str
    hash: str

    @staticmethod
    def from_sample_sheet(sample_sheet: SampleSheet) -> 'SampleSheetLookup':
        return SampleSheetLookup(sequencing_run=sample_sheet.sequencing_run.name, hash=sample_sheet.hash)


@dataclass_json
@dataclass
class Aligner:
    name: str
    version: str


@dataclass_json
@dataclass
class VariantCaller:
    name: str
    version: str
    run_params: Optional[str] = None


@dataclass_json
@dataclass
class SampleSheetCombinedVCFFile:
    path: str
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    variant_caller: VariantCaller


@dataclass_json
@dataclass
class BamFile:
    path: str
    aligner: Optional[Aligner] = field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class VCFFile:
    path: str
    variant_caller: Optional[VariantCaller] = field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class SequencingFile:
    sample_name: str
    fastq_r1: str
    fastq_r2: str
    bam_file: BamFile
    vcf_file: VCFFile


@dataclass_json
@dataclass
class SequencingSampleLookup:
    """ Only used as arguments to find existing sequencing sample on server - not enough details to create one """
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    sample_name: str


@dataclass_json
@dataclass
class QC:
    """ QC information needs to be matched against a particular BAM/VCF file (as there could be multiple)
        We use this to match gene lists, exec stats and coverage below
     """
    sequencing_sample_lookup: SequencingSampleLookup = field(metadata=config(field_name="sequencing_sample"))
    bam_file: BamFile
    vcf_file: VCFFile


@dataclass_json
@dataclass
class QCGeneList:
    path: str
    qc: QC
    gene_list: List[str]


@dataclass_json
@dataclass
class QCExecStats:
    path: str
    qc: QC
    created: datetime
    modified: datetime
    hash: str
    is_valid: bool
    mean_coverage_across_genes: float
    mean_coverage_across_kit: float
    median_insert: float
    percent_read_enrichment: float
    percent_duplication: float
    uniformity_of_coverage: float
    deduplicated_reads: Optional[int] = None
    indels_dbsnp_percent: Optional[float] = None
    number_indels: Optional[int] = None
    number_snps: Optional[int] = None
    percent_10x_goi: Optional[float] = None
    percent_20x_goi: Optional[float] = None
    percent_20x_kit: Optional[float] = None
    percent_100x_goi: Optional[float] = None
    percent_100x_kit: Optional[float] = None
    percent_250x_goi: Optional[float] = None
    percent_250x_kit: Optional[float] = None
    percent_500x_goi: Optional[float] = None
    percent_500x_kit: Optional[float] = None
    percent_error_rate: Optional[float] = None
    percent_map_to_diff_chr: Optional[float] = None
    percent_reads: Optional[float] = None
    percent_softclip: Optional[float] = None
    reads: Optional[int] = None
    sample_id_lod: Optional[float] = None
    sex_match: Optional[str] = None
    snp_dbsnp_percent: Optional[float] = None
    ts_to_tv_ratio: Optional[float] = None


@dataclass_json
@dataclass
class QCGeneCoverage:
    """ We send this up to associate coverage file path with sequencing sample
        Then, later we upload the coverage file - and match via path """
    qc: QC
    path: str
