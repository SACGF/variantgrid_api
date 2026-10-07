"""AlignmentFile / SequencingFile.alignment_files / QC.alignment_files (SACGF/variantgrid_api#27), and that code
written against the deprecated BamFile / bam_file still sends what it did before"""
import json

import pytest
import responses

from variantgrid_api.api_client import VariantGridAPI, UnsupportedFeaturePolicy, UnsupportedFeatureError
from variantgrid_api.data_models import (
    Aligner, AlignmentFile, AlignmentFileType, BamFile, QC, QCGeneList, SequencingFile, SequencingSampleLookup,
    SampleSheetLookup, ServerCapabilities, ServerFeature, SingleSampleVCF, VariantCaller,
)

ALIGNER = Aligner(name="BWA", version="0.7.18")
GATK = VariantCaller(name="GATK", version="4.1.9.0")
LOOKUP = SampleSheetLookup(sequencing_run="RUN_1", hash="abc")
SSL = SequencingSampleLookup(sample_sheet_lookup=LOOKUP, sample_name="s1")

# What 1.8.0 sent for a sample with FastQs, a BAM and a VCF - the shape older clients and servers rely on
OLD_RECORD = {
    "sample_name": "s1",
    "bam_file": {"path": "/data/s1.bam", "aligner": {"name": "BWA", "version": "0.7.18"}},
    "unaligned_reads": {"fastq_r1": {"path": "/data/s1_R1.fastq.gz"}, "fastq_r2": {"path": "/data/s1_R2.fastq.gz"}},
    "vcf_file": {"path": "/data/s1.vcf.gz", "variant_caller": {"name": "GATK", "version": "4.1.9.0", "run_params": None}},
}
OLD_QC_GENE_LIST = {
    "path": "/data/s1.goi.txt",
    "qc": {
        "sequencing_sample": {"sample_sheet": {"sequencing_run": "RUN_1", "hash": "abc"}, "sample_name": "s1"},
        "bam_file": {"path": "/data/s1.bam"},
        "vcf_file": {"path": "/data/s1.vcf.gz"},
    },
    "gene_list": ["BRCA1"],
}

SEQUENCING_FILES_URL = "https://example.org/seqauto/api/v1/sequencing_files/bulk_create"
QC_GENE_LIST_URL = "https://example.org/seqauto/api/v1/qc_gene_list/"
CAPABILITIES_URL = "https://example.org/api/v1/capabilities"


def _old_style_sequencing_file():
    with pytest.warns(DeprecationWarning) as record:
        sf = SequencingFile("s1", BamFile("/data/s1.bam", ALIGNER),
                            fastq_r1="/data/s1_R1.fastq.gz", fastq_r2="/data/s1_R2.fastq.gz",
                            vcf_files=[SingleSampleVCF("/data/s1.vcf.gz", GATK)])
    messages = {str(w.message) for w in record}
    assert "BamFile is deprecated; use AlignmentFile instead." in messages
    assert "SequencingFile.bam_file is deprecated; use alignment_files instead." in messages
    return sf


def _sequencing_file(*alignment_files):
    return SequencingFile("s1", alignment_files=list(alignment_files),
                          fastq_r1="/data/s1_R1.fastq.gz", fastq_r2="/data/s1_R2.fastq.gz",
                          vcf_files=[SingleSampleVCF("/data/s1.vcf.gz", GATK)])


def _api(features=None, **kwargs):
    api = VariantGridAPI("https://example.org", "TKN", **kwargs)
    if features is not None:
        api._capabilities = ServerCapabilities(version="test", features=frozenset(features))
    return api


def _post_sequencing_data(api, sequencing_files):
    responses.add(responses.POST, SEQUENCING_FILES_URL, json={"created": 1}, status=200)
    api.create_sequencing_data(LOOKUP, sequencing_files)
    return json.loads(responses.calls[-1].request.body)["records"]


def _post_qc_gene_list(api, qc):
    responses.add(responses.POST, QC_GENE_LIST_URL, json={}, status=200)
    api.create_qc_gene_list(QCGeneList(path="/data/s1.goi.txt", qc=qc, gene_list=["BRCA1"]))
    return json.loads(responses.calls[-1].request.body)


# ------------------------------------------------------------------ #
# Deprecated usage keeps working                                       #
# ------------------------------------------------------------------ #

def test_bam_file_is_deprecated_alias():
    with pytest.warns(DeprecationWarning, match="BamFile is deprecated"):
        bam_file = BamFile(path="/data/s1.bam", aligner=ALIGNER)
    assert isinstance(bam_file, AlignmentFile)
    assert bam_file.to_dict() == {"path": "/data/s1.bam", "aligner": {"name": "BWA", "version": "0.7.18"}}
    assert bam_file.to_dict() == AlignmentFile(path="/data/s1.bam", aligner=ALIGNER).to_dict()


@responses.activate
def test_old_style_sequencing_file_sends_1_8_0_json_without_probe():
    """ bam_file alone is sent exactly as before, and never asks the server for its capabilities """
    sf = _old_style_sequencing_file()
    assert sf.bam_file.path == "/data/s1.bam"
    assert sf.get_alignment_files() == [sf.bam_file]
    assert _post_sequencing_data(_api(), [sf]) == [OLD_RECORD]
    assert len(responses.calls) == 1


@responses.activate
def test_old_style_qc_sends_1_8_0_json_without_probe():
    with pytest.warns(DeprecationWarning, match="QC.bam_file is deprecated"):
        qc = QC(SSL, AlignmentFile("/data/s1.bam"), SingleSampleVCF("/data/s1.vcf.gz"))  # positional, as before
    assert qc.get_alignment_files() == [qc.bam_file]
    assert _post_qc_gene_list(_api(), qc) == OLD_QC_GENE_LIST
    assert len(responses.calls) == 1


@responses.activate
def test_qc_bam_file_merged_into_alignment_files():
    with pytest.warns(DeprecationWarning, match="QC.bam_file is deprecated"):
        qc = QC(SSL, bam_file=AlignmentFile("/data/s1.bam"), vcf_file=SingleSampleVCF("/data/s1.vcf.gz"),
                alignment_files=[AlignmentFile("/data/s1.cram")])
    assert [af.path for af in qc.get_alignment_files()] == ["/data/s1.bam", "/data/s1.cram"]
    body = _post_qc_gene_list(_api(features=[ServerFeature.ALIGNMENT_FILES]), qc)
    assert "bam_file" not in body["qc"]
    assert body["qc"]["alignment_files"] == [{"path": "/data/s1.bam"}, {"path": "/data/s1.cram"}]


# ------------------------------------------------------------------ #
# alignment_files                                                      #
# ------------------------------------------------------------------ #

def test_alignment_file_type_inferred_from_path_and_not_sent():
    assert AlignmentFile("/data/S1.CRAM").get_file_type() == AlignmentFileType.CRAM
    assert AlignmentFile("/data/s1.bam").get_file_type() == AlignmentFileType.BAM
    assert AlignmentFile("/data/s1.cram", file_type=AlignmentFileType.BAM).get_file_type() == AlignmentFileType.BAM
    assert "file_type" not in AlignmentFile("/data/s1.cram").to_dict()


@responses.activate
def test_alignment_files_sent_to_server_with_feature():
    sf = _sequencing_file(AlignmentFile("/data/s1.bam", ALIGNER),
                          AlignmentFile("/data/s1.cram", ALIGNER, file_type=AlignmentFileType.CRAM))
    records = _post_sequencing_data(_api(features=[ServerFeature.ALIGNMENT_FILES]), [sf])
    assert len(records) == 1
    record = records[0]
    assert "bam_file" not in record
    assert [af["path"] for af in record["alignment_files"]] == ["/data/s1.bam", "/data/s1.cram"]
    assert record["alignment_files"][1]["file_type"] == "cram"


@responses.activate
def test_deprecated_bam_file_merged_into_alignment_files():
    with pytest.warns(DeprecationWarning):
        sf = SequencingFile("s1", bam_file=AlignmentFile("/data/s1.bam"),
                            alignment_files=[AlignmentFile("/data/s1.recal.bam")],
                            vcf_files=[SingleSampleVCF("/data/s1.vcf.gz", GATK)])
    record = _post_sequencing_data(_api(features=[ServerFeature.ALIGNMENT_FILES]), [sf])[0]
    assert "bam_file" not in record
    assert [af["path"] for af in record["alignment_files"]] == ["/data/s1.bam", "/data/s1.recal.bam"]


@responses.activate
def test_legacy_server_one_record_per_alignment_file_as_bam_file():
    """ An older server gets what 1.8.0 sent - a single BAM is the identical record """
    responses.add(responses.GET, CAPABILITIES_URL, status=404)
    api = _api()
    assert _post_sequencing_data(api, [_sequencing_file(AlignmentFile("/data/s1.bam", ALIGNER))]) == [OLD_RECORD]

    sf = _sequencing_file(AlignmentFile("/data/s1.bam", ALIGNER, file_type=AlignmentFileType.BAM),
                          AlignmentFile("/data/s1.recal.bam", ALIGNER))
    records = _post_sequencing_data(api, [sf])
    assert [r["bam_file"]["path"] for r in records] == ["/data/s1.bam", "/data/s1.recal.bam"]
    assert records[0] == OLD_RECORD  # file_type isn't sent as bam_file
    assert all("alignment_files" not in r for r in records)


@responses.activate
def test_server_without_alignment_files_sends_cram_as_bam_file():
    sf = _sequencing_file(AlignmentFile("/data/s1.cram"))
    records = _post_sequencing_data(_api(features=[ServerFeature.CRAM_ALIGNMENT_FILES]), [sf])
    assert [r["bam_file"] for r in records] == [{"path": "/data/s1.cram"}]


def test_cram_to_server_without_cram_raises():
    sf = _sequencing_file(AlignmentFile("/data/s1.bam"), AlignmentFile("/data/s1.cram"))
    with pytest.raises(UnsupportedFeatureError, match="cram_alignment_files"):
        _api(features=[]).create_sequencing_data(LOOKUP, [sf])


@responses.activate
def test_cram_to_server_without_cram_skipped_under_skip():
    sf = _sequencing_file(AlignmentFile("/data/s1.bam"), AlignmentFile("/data/s1.cram"))
    api = _api(features=[], unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)
    records = _post_sequencing_data(api, [sf])
    assert [r["bam_file"]["path"] for r in records] == ["/data/s1.bam"]


@responses.activate
def test_qc_alignment_files_sent_to_server_with_feature():
    qc = QC(SSL, vcf_file=SingleSampleVCF("/data/s1.vcf.gz"),
            alignment_files=[AlignmentFile("/data/s1.bam"), AlignmentFile("/data/s1.cram")])
    body = _post_qc_gene_list(_api(features=[ServerFeature.ALIGNMENT_FILES]), qc)
    assert body["qc"]["alignment_files"] == [{"path": "/data/s1.bam"}, {"path": "/data/s1.cram"}]
    assert "bam_file" not in body["qc"]


@responses.activate
def test_qc_first_alignment_file_sent_as_bam_file_to_legacy_server():
    """ An older server finds the QC by sample, VCF path and bam_file path - the first alignment file """
    qc = QC(SSL, vcf_file=SingleSampleVCF("/data/s1.vcf.gz"),
            alignment_files=[AlignmentFile("/data/s1.bam", file_type=AlignmentFileType.BAM),
                             AlignmentFile("/data/s1.cram")])
    assert _post_qc_gene_list(_api(features=[]), qc) == OLD_QC_GENE_LIST
