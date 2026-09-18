import re
import warnings
from datetime import date, datetime
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, List, Union, FrozenSet, ClassVar

from dataclasses_json import dataclass_json, config


@dataclass_json
@dataclass
class EnrichmentKit:
    name: str
    version: int


@dataclass_json
@dataclass
class Manufacturer:
    name: str

@dataclass_json
@dataclass
class SequencerModel:
    model: str
    manufacturer: Manufacturer
    data_naming_convention: str # 'M' or 'H' for (HiSeq only)

    # (prefix, model_name, data_naming_convention)
    _ILLUMINA_PREFIXES = [
        ("LH", "NovaSeq X",    "M"),
        ("NB", "NextSeq",      "M"),
        ("NS", "NextSeq",      "M"),
        ("A",  "NovaSeq 6000", "M"),
        ("H",  "HiSeq",        "H"),
        ("M",  "MiSeq",        "M"),
    ]

    @classmethod
    def from_sequencer_name(cls, sequencer_name: str) -> 'SequencerModel':
        """Infer SequencerModel from an Illumina sequencer instrument name.

        Illumina instrument names begin with a letter code that identifies the
        platform (e.g. ``M02027`` → MiSeq, ``A01234`` → NovaSeq 6000).
        Returns an Unknown model when no prefix matches.
        """
        illumina = Manufacturer(name="Illumina")
        for prefix, model_name, convention in cls._ILLUMINA_PREFIXES:
            if sequencer_name.startswith(prefix):
                return cls(model=model_name, manufacturer=illumina, data_naming_convention=convention)
        return cls(model="Unknown", manufacturer=illumina, data_naming_convention="U")


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
class SequencingSampleLookup:
    """ Only used as arguments to find existing sequencing sample on server - not enough details to create one """
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    sample_name: str


@dataclass_json
@dataclass
class JointCalledVCF:
    """A joint-called multi-sample VCF (mirrors the server ``JointCalledVCF`` model). """
    path: str
    sample_sheet_lookup: SampleSheetLookup = field(metadata=config(field_name="sample_sheet"))
    variant_caller: VariantCaller
    # Set for joint calls that draw samples from more than one sequencing run, eg a family trio.
    # sample_sheet_lookup stays the owning run - the one the path sits under.
    sequencing_samples: Optional[List[SequencingSampleLookup]] = \
        field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class SampleSheetCombinedVCFFile(JointCalledVCF):
    """Deprecated alias for :class:`JointCalledVCF` - use that instead. """
    def __post_init__(self):
        warnings.warn(
            "SampleSheetCombinedVCFFile is deprecated; use JointCalledVCF instead.",
            DeprecationWarning,
            stacklevel=2,
        )


@dataclass_json
@dataclass
class BamFile:
    path: str
    aligner: Optional[Aligner] = field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class SingleSampleVCF:
    """A per-sample VCF - one BAM in, one VCF out (mirrors the server ``SingleSampleVCF`` model). """
    path: str
    variant_caller: Optional[VariantCaller] = field(default=None, metadata=config(exclude=lambda x: x is None))


@dataclass_json
@dataclass
class VCFFile(SingleSampleVCF):
    """Deprecated alias for :class:`SingleSampleVCF` - use that instead. """
    def __post_init__(self):
        warnings.warn(
            "VCFFile is deprecated; use SingleSampleVCF instead.",
            DeprecationWarning,
            stacklevel=2,
        )


@dataclass_json
@dataclass
class SequencingFile:
    """ FastQs are optional - BAM-first runs (sequencer emits BAM, or FastQs not kept) send just BAM + VCF

        vcf_files: the VCFs called off this BAM, one per variant caller, eg DRAGEN TSO 500's small variant VCF and
        its gene-level CNV VCF. The server keeps one VCF per BAM and caller, so a second with the same caller
        would replace the first - create_sequencing_data() raises instead.

        vcf_file is deprecated - use vcf_files. It still works (set it and it is sent, read it back as before),
        and get_vcf_files() gives both """
    sample_name: str
    bam_file: BamFile
    vcf_file: Optional[SingleSampleVCF] = field(default=None, metadata=config(exclude=lambda x: x is None))
    fastq_r1: Optional[str] = field(default=None, metadata=config(exclude=lambda x: x is None))
    fastq_r2: Optional[str] = field(default=None, metadata=config(exclude=lambda x: x is None))
    vcf_files: Optional[List[SingleSampleVCF]] = field(default=None, metadata=config(exclude=lambda x: x is None))

    def __post_init__(self):
        if self.vcf_file is not None:
            warnings.warn(
                "SequencingFile.vcf_file is deprecated; use vcf_files instead.",
                DeprecationWarning,
                stacklevel=3,
            )

    def get_vcf_files(self) -> List[SingleSampleVCF]:
        """ vcf_file (deprecated) then vcf_files """
        vcf_files = [self.vcf_file] if self.vcf_file is not None else []
        return vcf_files + list(self.vcf_files or [])


@dataclass_json
@dataclass
class QC:
    """ QC information needs to be matched against a particular BAM/VCF file (as there could be multiple)
        We use this to match gene lists, exec stats and coverage below
     """
    sequencing_sample_lookup: SequencingSampleLookup = field(metadata=config(field_name="sequencing_sample"))
    bam_file: BamFile
    vcf_file: SingleSampleVCF


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


############################################################
## Patient -> Specimen -> Extraction (SACGF/variantgrid#1707)

class Sex(str, Enum):
    UNKNOWN = "U"
    MALE = "M"
    FEMALE = "F"


class TissueStatus(str, Enum):
    """ The role the material plays in the test. Blood is the reference in a solid tumour workup
        but the tumour in leukaemia, so this belongs to the specimen, not the tissue """
    REFERENCE = "R"  # reference / unaffected
    AFFECTED = "A"   # affected / lesional
    UNKNOWN = "U"


class NucleicAcid(str, Enum):
    DNA = "D"
    RNA = "R"


class SpecimenMeasureType(str, Enum):
    TMB = "T"              # Tumour mutational burden
    MSI = "M"              # Microsatellite instability
    GIS = "G"              # Genomic instability score
    TUMOUR_FRACTION = "F"
    PLOIDY = "P"


def _exclude_none():
    return config(exclude=lambda x: x is None)


@dataclass_json
@dataclass
class ExternalPK:
    """ The key an external system (eg a LIMS) knows a record by. The server's ExternalPK is unique on all
        three, so all three are required """
    code: str
    external_type: str      # which identifier scheme the code belongs to, eg 'HelixID'
    external_manager: str   # the system that owns the record, eg 'HELIX'


@dataclass_json
@dataclass
class ExternalReference:
    """ Names an existing Patient / Specimen / Extraction: by its local reference (reference_id, or
        patient_code for a patient), by an ExternalPK (code + external_type, optionally narrowed to
        one external_manager), or by both.

        Wherever a reference is taken a bare string works too, and means the local reference. """
    reference_id: Optional[str] = field(default=None, metadata=_exclude_none())
    code: Optional[str] = field(default=None, metadata=_exclude_none())
    external_type: Optional[str] = field(default=None, metadata=_exclude_none())
    external_manager: Optional[str] = field(default=None, metadata=_exclude_none())

    def __post_init__(self):
        if not (self.reference_id or self.code):
            raise ValueError("ExternalReference must supply 'reference_id' or 'code'")
        if bool(self.code) != bool(self.external_type):
            # A code alone names nothing - the server would 400
            raise ValueError("ExternalReference 'code' and 'external_type' must be supplied together")

    def to_json_value(self) -> Union[str, dict]:
        """ A bare string when only reference_id is set, otherwise an object """
        data = self.to_dict()
        if list(data) == ["reference_id"]:
            return self.reference_id
        return data

    @staticmethod
    def _from_local_and_external(local_reference: Optional[str],
                                 external_pk: Optional[ExternalPK]) -> 'ExternalReference':
        kwargs = {"reference_id": local_reference}
        if external_pk:
            kwargs.update(code=external_pk.code, external_type=external_pk.external_type,
                          external_manager=external_pk.external_manager)
        return ExternalReference(**kwargs)

    @staticmethod
    def from_patient(patient: 'Patient') -> 'ExternalReference':
        return ExternalReference._from_local_and_external(patient.patient_code, patient.external_pk)

    @staticmethod
    def from_specimen(specimen: 'Specimen') -> 'ExternalReference':
        return ExternalReference._from_local_and_external(specimen.reference_id, specimen.external_pk)

    @staticmethod
    def from_extraction(extraction: 'Extraction') -> 'ExternalReference':
        return ExternalReference._from_local_and_external(extraction.reference_id, extraction.external_pk)


ReferenceLike = Union[str, ExternalReference]


def reference_json(reference: Optional[ReferenceLike]) -> Union[str, dict, None]:
    """ JSON for a reference given as a bare string or an ExternalReference """
    if isinstance(reference, ExternalReference):
        return reference.to_json_value()
    return reference


def _reference_field(**kwargs):
    return field(metadata=config(encoder=reference_json, exclude=lambda x: x is None), **kwargs)


@dataclass_json
@dataclass
class Patient:
    """ Everything is optional, but send patient_code and/or external_pk - they are what a re-post
        matches on, so without either every post creates a new patient """
    patient_code: Optional[str] = field(default=None, metadata=_exclude_none())
    family_code: Optional[str] = field(default=None, metadata=_exclude_none())
    first_name: Optional[str] = field(default=None, metadata=_exclude_none())
    last_name: Optional[str] = field(default=None, metadata=_exclude_none())
    date_of_birth: Optional[date] = field(default=None, metadata=_exclude_none())
    date_of_death: Optional[date] = field(default=None, metadata=_exclude_none())
    sex: Optional[Sex] = field(default=None, metadata=_exclude_none())
    affected: Optional[bool] = field(default=None, metadata=_exclude_none())
    external_pk: Optional[ExternalPK] = field(default=None, metadata=_exclude_none())


@dataclass_json
@dataclass
class Specimen:
    """ reference_id is unique per patient, not globally """
    patient: ReferenceLike = _reference_field()
    reference_id: str
    description: Optional[str] = field(default=None, metadata=_exclude_none())
    collected_by: Optional[str] = field(default=None, metadata=_exclude_none())
    collection_date: Optional[datetime] = field(default=None, metadata=_exclude_none())
    received_date: Optional[datetime] = field(default=None, metadata=_exclude_none())
    tissue_status: Optional[TissueStatus] = field(default=None, metadata=_exclude_none())
    external_pk: Optional[ExternalPK] = field(default=None, metadata=_exclude_none())


@dataclass_json
@dataclass
class Extraction:
    """ One nucleic acid extraction off a specimen - a TSO 500 specimen has two, the DNA and RNA arms """
    specimen: ReferenceLike = _reference_field()
    reference_id: Optional[str] = field(default=None, metadata=_exclude_none())
    nucleic_acid_source: Optional[NucleicAcid] = field(default=None, metadata=_exclude_none())
    extraction_date: Optional[datetime] = field(default=None, metadata=_exclude_none())
    external_pk: Optional[ExternalPK] = field(default=None, metadata=_exclude_none())


@dataclass_json
@dataclass
class SpecimenMeasure:
    """ A scalar measured on the specimen rather than on any one variant - TMB, MSI, GIS etc (SACGF/variantgrid#1559)

        Transcribe these from vendor output rather than computing them, and put the raw block they came
        from in source_payload. Send both the score (value) and the lab's call, with the threshold that
        turned one into the other - that threshold is lab policy, not vendor output.

        The specimen is passed to create_specimen_measure(s) rather than held here. A re-post for the
        same specimen and measure_type replaces the previous value. """
    measure_type: SpecimenMeasureType
    value: Optional[float] = field(default=None, metadata=_exclude_none())
    unit: Optional[str] = field(default=None, metadata=_exclude_none())             # eg 'mut/Mb', '%'
    call: Optional[str] = field(default=None, metadata=_exclude_none())             # eg 'High', 'Stable'
    threshold: Optional[str] = field(default=None, metadata=_exclude_none())
    threshold_source: Optional[str] = field(default=None, metadata=_exclude_none())  # whose policy set it
    method: Optional[str] = field(default=None, metadata=_exclude_none())           # tool and version
    source_payload: Optional[dict] = field(default=None, metadata=_exclude_none())
    measured_date: Optional[datetime] = field(default=None, metadata=_exclude_none())
    extraction: Optional[ReferenceLike] = _reference_field(default=None)  # the arm that produced it


############################################################
## Server capabilities (SACGF/variantgrid_api#22)

class _ServerName(str, Enum):
    """ A name from the server's vocabulary. Subclasses str so plain strings and these are interchangeable,
        and formats as its value so messages read 'patients' rather than 'ServerFeature.PATIENTS' """

    def __str__(self):
        return self.value

    def __format__(self, format_spec):
        return format(self.value, format_spec)


class ServerFeature(_ServerName):
    """ Features a server reports in capabilities (API_FEATURES in the variantgrid repo's
        variantgrid/views_rest.py). Names are never removed - an older server just doesn't list a newer one """
    PATIENTS = "patients"
    SPECIMEN_MEASURES = "specimen_measures"
    LINK_EXTRACTION = "link_extraction"
    UPLOAD_STATUS = "upload_status"
    JOINT_CALLED_VCF_CROSS_RUN = "joint_called_vcf_cross_run"
    UPLOAD_METADATA = "upload_metadata"


class UploadFileType(_ServerName):
    """ File types a server imports from an upload - the server's upload.models.UploadedFileTypes names in
        lower case, less the internal ones it drives itself. A server only reports those it has an importer for """
    BED = "bed"
    DRAGEN_TSO500_ALL_FUSIONS = "dragen_tso500_all_fusions"
    DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT = "dragen_tso500_combined_variant_output"
    GENE_COVERAGE = "gene_coverage"
    GENE_LIST = "gene_list"
    GENE_LEVEL_CNV_VCF = "gene_level_cnv_vcf"
    GENE_LEVEL_INSERT_VARIANTS_ONLY = "gene_level_insert_variants_only"
    PATIENT_RECORDS = "patient_records"
    PED = "ped"
    VCF = "vcf"
    VCF_INSERT_VARIANTS_ONLY = "vcf_insert_variants_only"


@dataclass(frozen=True)
class ServerCapabilities:
    """ What a server accepts, from GET api/v1/capabilities (SACGF/variantgrid_sapath#443).

        A server without that endpoint (404) is LEGACY: version 'legacy', no features, no upload file types """
    LEGACY: ClassVar["ServerCapabilities"]

    version: str
    git_hash: Optional[str] = None
    features: FrozenSet[str] = frozenset()
    upload_file_types: FrozenSet[str] = frozenset()

    @staticmethod
    def from_json(data: dict) -> "ServerCapabilities":
        return ServerCapabilities(version=data["version"],
                                  git_hash=data.get("git_hash"),
                                  features=frozenset(data.get("features") or []),
                                  upload_file_types=frozenset(data.get("upload_file_types") or []))


ServerCapabilities.LEGACY = ServerCapabilities(version="legacy")
