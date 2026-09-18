import datetime
import json
import logging
import re
import time
import urllib
import warnings
from enum import Enum
from pathlib import Path
from typing import List, Optional, Callable, Union

import requests

from variantgrid_api.data_models import EnrichmentKit, SequencingRun, SampleSheet, JointCalledVCF, \
    SampleSheetLookup, SequencingFile, QCGeneList, QCExecStats, QCGeneCoverage, SequencerModel, Sequencer, \
    SequencingSampleLookup, Patient, Specimen, Extraction, SpecimenMeasure, ExternalReference, ReferenceLike, \
    reference_json, ServerCapabilities, ServerFeature, UploadFileType


_UNSET = object()


class DateTimeEncoder(json.JSONEncoder):
    def default(self,o):
        if isinstance(o,(datetime.date, datetime.datetime)):
            return o.isoformat()
        return super().default(o)


class EmptyInputPolicy(Enum):
    IGNORE = "ignore"
    WARN = "warn"
    ERROR = "error"


class UnsupportedFeaturePolicy(Enum):
    """ What a gated call does when the server lacks the feature (see VariantGridAPI.capabilities) """
    SKIP = "skip"     # log a warning, return None
    ERROR = "error"   # raise UnsupportedFeatureError


class UnsupportedFeatureError(Exception):
    """ Raised (under UnsupportedFeaturePolicy.ERROR) when a call needs a feature the server doesn't have """

    def __init__(self, message, capabilities: Optional[ServerCapabilities] = None):
        super().__init__(message)
        self.capabilities = capabilities


class AnnotationError(Exception):
    """Raised when the server reports an error while importing/annotating an uploaded file."""

    def __init__(self, message, status: Optional[dict] = None):
        super().__init__(message)
        self.status = status


class VariantGridAPI:
    def __init__(self, server, api_token,
                 empty_input_policy=EmptyInputPolicy.ERROR,
                 logger: Optional[logging.Logger] = None,
                 log_request=False, log_response=False,
                 unsupported_feature_policy=UnsupportedFeaturePolicy.ERROR
    ):
        self.server = server
        self.headers = {"Authorization": f"Token {api_token}"}
        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger
        self.validation_handler = self._get_validation_handler(empty_input_policy, logger)
        self.log_request = log_request
        self.log_response = log_response
        self.unsupported_feature_policy = unsupported_feature_policy
        self._capabilities: Optional[ServerCapabilities] = None

    def _get_url(self, url):
        return urllib.parse.urljoin(self.server, url)

    def _post(self, path, json_data):
        url = self._get_url(path)
        json_string = json.dumps(json_data, cls=DateTimeEncoder)
        if self.log_request:
            self.logger.info("POST to '%s', JSON: %s", url, json_string)
        response = requests.post(url,
                                 headers={**self.headers, "Content-Type": "application/json"},
                                 data=json_string)

        extra_error_message = f"{url=}, {json_string=}"
        return self._handle_json_response(response, extra_error_message)

    def _handle_json_response(self, response, extra_error_message: Optional[str] = None):
        try:
            json_response = response.json()
            if self.log_response:
                self.logger.info("Response from '%s', JSON: %s", response.url, json_response)
        except Exception as e:
            json_response = f"Couldn't convert JSON: {e}"
        if not response.ok:
            if extra_error_message:
                self.logger.error(extra_error_message)
            self.logger.error("Response: %s", json_response)
            response.raise_for_status()
        return json_response

    @staticmethod
    def _get_validation_handler(empty_input_policy: EmptyInputPolicy, logger: logging.Logger) -> Callable:
        def _ignore(msg):
            pass

        def _warn(msg):
            logger.warning(msg)

        def _error(msg):
            raise ValueError(msg)

        return {
            EmptyInputPolicy.IGNORE: _ignore,
            EmptyInputPolicy.WARN: _warn,
            EmptyInputPolicy.ERROR: _error,
        }[empty_input_policy]

    def _validate_string(self, info: str, s: Optional[str]) -> None:
        if s is None:
            self.validation_handler(f"{info}: is None")
        elif s == "":
            self.validation_handler(f"{info}: empty string")

    def _validate_object(self, info: str, obj):
        if obj is None:
            self.validation_handler(f"{info}: is None")

    def _validate_reference(self, info: str, reference: Optional[ReferenceLike]):
        if isinstance(reference, str):
            self._validate_string(info, reference)
        else:
            self._validate_object(info, reference)

    def _validate_list(self, info: str, list_obj: List):
        if not list_obj:
            self.validation_handler(f"{info}: empty list")

    ##########################################################
    ## Server capabilities (SACGF/variantgrid_sapath#443)
    ## One client talks to servers of different ages. Calls needing a newer server are gated on the
    ## features it reports, and handled by unsupported_feature_policy when it lacks them

    @property
    def capabilities(self) -> ServerCapabilities:
        """ Fetched on first use then cached, so constructing the client and ungated calls do no extra
            request. A server without the endpoint is ServerCapabilities.LEGACY - it answers 404, or a redirect
            to its login page when its login middleware doesn't exempt /api/ (so redirects aren't followed) """
        if self._capabilities is None:
            url = self._get_url("api/v1/capabilities")
            response = requests.get(url, headers=self.headers, allow_redirects=False)
            if response.status_code == 404 or response.is_redirect:
                self.logger.info("Server has no capabilities endpoint (HTTP %s) - treating as legacy",
                                 response.status_code)
                self._capabilities = ServerCapabilities.LEGACY
            else:
                data = self._handle_json_response(response, f"{url=}")
                self._capabilities = ServerCapabilities.from_json(data)
        return self._capabilities

    def supports(self, feature: Union[ServerFeature, str]) -> bool:
        return feature in self.capabilities.features

    def accepts_upload(self, file_type: Union[UploadFileType, str]) -> bool:
        """ file_type is an UploadFileType, or the server's UploadedFileTypes name in lower case """
        return file_type in self.capabilities.upload_file_types

    def _unsupported(self, message: str) -> bool:
        """ Applies unsupported_feature_policy - returns False (caller skips) or raises """
        capabilities = self.capabilities
        message = f"{message} (server version '{capabilities.version}')"
        if self.unsupported_feature_policy == UnsupportedFeaturePolicy.SKIP:
            self.logger.warning("Skipping: %s", message)
            return False
        raise UnsupportedFeatureError(message, capabilities)

    def _require(self, feature: Union[ServerFeature, str]) -> bool:
        """ True if the server supports feature, otherwise applies unsupported_feature_policy """
        return self.supports(feature) or self._unsupported(f"server doesn't support feature '{feature}'")

    def _require_upload(self, file_type: Union[UploadFileType, str]) -> bool:
        return self.accepts_upload(file_type) or self._unsupported(f"server doesn't accept upload file type '{file_type}'")

    def create_experiment(self, experiment: str):
        self._validate_string("experiment", experiment)
        json_data = {
            "name": experiment
        }
        return self._post("seqauto/api/v1/experiment/", json_data)

    def create_enrichment_kit(self, enrichment_kit: EnrichmentKit):
        self._validate_object("enrichment_kit", enrichment_kit)
        return self._post("seqauto/api/v1/enrichment_kit/",
                          enrichment_kit.to_dict())


    def create_sequencer_model(self, sequencer_model: SequencerModel):
        self._validate_object("sequencer_model", sequencer_model)
        return self._post("seqauto/api/v1/sequencer_model/",
                          sequencer_model.to_dict())

    def create_sequencer(self, sequencer: Sequencer):
        self._validate_object("sequencer", sequencer)
        return self._post("seqauto/api/v1/sequencer/",
                          sequencer.to_dict())

    def create_sequencing_run(self, sequencing_run: SequencingRun):
        self._validate_object("sequencing_run", sequencing_run)
        return self._post("seqauto/api/v1/sequencing_run/",
                          sequencing_run.to_dict())

    def create_sample_sheet(self, sample_sheet: SampleSheet):
        self._validate_object("sample_sheet", sample_sheet)
        json_data = sample_sheet.to_dict()
        # We don't want all sequencing_run just the name
        sequencing_run = json_data.pop("sequencing_run")
        json_data["sequencing_run"] = sequencing_run["name"]
        return self._post("seqauto/api/v1/sample_sheet/",
                          json_data)

    def create_joint_called_vcf(self, joint_called_vcf: JointCalledVCF):
        self._validate_object("joint_called_vcf", joint_called_vcf)
        json_data = joint_called_vcf.to_dict()
        return self._post("seqauto/api/v1/joint_called_vcf/",
                          json_data)

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file: JointCalledVCF):
        """Deprecated alias for :meth:`create_joint_called_vcf` - use that instead. """
        warnings.warn(
            "create_sample_sheet_combined_vcf_file is deprecated; use create_joint_called_vcf instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.create_joint_called_vcf(sample_sheet_combined_vcf_file)

    def create_sequencing_data(self, sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]):
        self._validate_object("sample_sheet_lookup", sample_sheet_lookup)
        self._validate_list("sequencing_files", sequencing_files)
        records = []
        for sf in sequencing_files:
            # The server requires both paths - catch it here, naming the record, rather than a 400 for the batch
            self._validate_string(f"SequencingFile '{sf.sample_name}' bam_file.path", sf.bam_file and sf.bam_file.path)
            vcf_files = sf.get_vcf_files()
            self._validate_list(f"SequencingFile '{sf.sample_name}' vcf_files", vcf_files)
            for i, vcf_file in enumerate(vcf_files):
                self._validate_string(f"SequencingFile '{sf.sample_name}' vcf_files[{i}].path",
                                      vcf_file and vcf_file.path)
            # The server keeps one VCF per BAM and caller - a repeated caller would silently replace a path
            callers = [f"{vc.name} {vc.version}" for vcf_file in vcf_files
                       if vcf_file and (vc := vcf_file.variant_caller)]
            if repeated := {c for c in callers if callers.count(c) > 1}:
                raise ValueError(f"SequencingFile '{sf.sample_name}' has more than one VCF from variant caller(s) "
                                 f"{', '.join(sorted(repeated))} - each VCF off a BAM needs its own caller")
            data = sf.to_dict()
            data.pop("vcf_file", None)
            data.pop("vcf_files", None)
            # put into hierarchial JSON DRF expects
            fastq_r1 = data.pop("fastq_r1", None)
            fastq_r2 = data.pop("fastq_r2", None)
            if fastq_r1:
                unaligned_reads = {"fastq_r1": {"path": fastq_r1}}
                if fastq_r2:
                    unaligned_reads["fastq_r2"] = {"path": fastq_r2}
                data["unaligned_reads"] = unaligned_reads
            elif fastq_r2:
                raise ValueError(f"SequencingFile '{sf.sample_name}' has fastq_r2 without fastq_r1")
            # No FastQs (BAM-first run) - server resolves the sample from sample_name.
            # The server takes one VCF per record, so each is a record sharing the BAM and FastQs
            for vcf_file in vcf_files:
                records.append({**data, "vcf_file": vcf_file.to_dict() if vcf_file else None})

        json_data = {
            "sample_sheet": sample_sheet_lookup.to_dict(),
            "records": records
        }
        return self._post("seqauto/api/v1/sequencing_files/bulk_create",
                          json_data)

    def create_qc_gene_list(self, qc_gene_list: QCGeneList):
        self._validate_object("qc_gene_list", qc_gene_list)
        json_data = qc_gene_list.to_dict()
        return self._post("seqauto/api/v1/qc_gene_list/",
                          json_data)


    def create_multiple_qc_gene_lists(self, qc_gene_lists: List[QCGeneList]):
        self._validate_list("qc_gene_lists", qc_gene_lists)
        json_data = {
            "records": [
                qcgl.to_dict() for qcgl in qc_gene_lists
            ]
        }
        return self._post("seqauto/api/v1/qc_gene_list/bulk_create",
                          json_data)

    def create_qc_exec_stats(self, qc_exec_stats: QCExecStats):
        self._validate_object("qc_exec_stats", qc_exec_stats)
        json_data = qc_exec_stats.to_dict()
        return self._post("seqauto/api/v1/qc_exec_summary/",
                          json_data)

    def create_multiple_qc_exec_stats(self, qc_exec_stats: List[QCExecStats]):
        self._validate_list("qc_exec_stats", qc_exec_stats)
        json_data = {
            "records": [
                qces.to_dict() for qces in qc_exec_stats
            ]
        }
        return self._post("seqauto/api/v1/qc_exec_summary/bulk_create",
                          json_data)

    def create_multiple_qc_gene_coverage(self, qc_gene_coverage_list: List[QCGeneCoverage]):
        self._validate_list("qc_gene_coverage_list", qc_gene_coverage_list)
        json_data = {
            "records": [
                qcgc.to_dict() for qcgc in qc_gene_coverage_list
            ]
        }
        return self._post("seqauto/api/v1/qc_gene_coverage/bulk_create",
                          json_data)

    ##########################################################
    ## Patient -> Specimen -> Extraction (SACGF/variantgrid#1707)
    ## Creates are upserts keyed on the identifiers sent, so re-posting returns the same rows

    def create_patient(self, patient: Patient):
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._validate_object("patient", patient)
        return self._post("patients/api/v1/patient/", patient.to_dict())

    def create_specimen(self, specimen: Specimen):
        """ The specimen's patient must already exist on the server, otherwise this is a 400 """
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._validate_object("specimen", specimen)
        return self._post("patients/api/v1/specimen/", specimen.to_dict())

    def create_extraction(self, extraction: Extraction):
        """ The extraction's specimen must already exist on the server, otherwise this is a 400 """
        if not self._require(ServerFeature.PATIENTS):
            return None
        self._validate_object("extraction", extraction)
        return self._post("patients/api/v1/extraction/", extraction.to_dict())

    def create_specimen_measure(self, specimen_reference: ReferenceLike, measure: SpecimenMeasure):
        """ An unknown specimen is a 400. Replaces any existing measure of the same type for the specimen """
        if not self._require(ServerFeature.SPECIMEN_MEASURES):
            return None
        self._validate_reference("specimen_reference", specimen_reference)
        self._validate_object("measure", measure)
        json_data = {"specimen": reference_json(specimen_reference), **measure.to_dict()}
        return self._post("patients/api/v1/specimen_measure/", json_data)

    def create_specimen_measures(self, specimen_reference: ReferenceLike, measures: List[SpecimenMeasure]):
        """ A run's measures (TMB, MSI, GIS etc) against one specimen in one call """
        if not self._require(ServerFeature.SPECIMEN_MEASURES):
            return None
        self._validate_reference("specimen_reference", specimen_reference)
        self._validate_list("measures", measures)
        json_data = {
            "specimen": reference_json(specimen_reference),
            "measures": [measure.to_dict() for measure in measures],
        }
        return self._post("patients/api/v1/specimen_measure/bulk_create", json_data)

    def link_sequencing_sample_extraction(self, sequencing_sample_lookup: SequencingSampleLookup,
                                          extraction_reference: ReferenceLike) -> Optional[dict]:
        """ Name the extraction a sequencing sample was made from. One call per sequencing sample is
            enough - the server carries it to every Sample made from that sample's VCFs, and on to the
            new rows if the sample sheet is re-sent.

            Returns {"sequencing_sample", "match_status", "match_error", "extraction"}. An unknown
            sequencing sample is a 400, but an extraction the server doesn't have yet is not an error:
            the response is a 202 with match_status 'Pending', and the link attaches itself once the
            extraction is created - there's no need to re-send. """
        if not self._require(ServerFeature.LINK_EXTRACTION):
            return None
        self._validate_object("sequencing_sample_lookup", sequencing_sample_lookup)
        self._validate_reference("extraction_reference", extraction_reference)
        json_data = {
            "sequencing_sample": sequencing_sample_lookup.to_dict(),
            "extraction": reference_json(extraction_reference),
        }
        return self._post("seqauto/api/v1/sequencing_sample/link_extraction", json_data)

    def upload_file(self, filename: str, path=_UNSET, metadata: Optional[dict] = None,
                    file_type: Optional[Union[UploadFileType, str]] = None):
        """ Upload a file via multipart POST to upload/api/v1/file_upload.

            Returns {"uploaded_file_id": <id>, "sha256_hash": <hash>, ...}; identify the upload by
            uploaded_file_id and/or sha256_hash (the server dedups on the content hash).

            path: SeqAuto-only server-side hint that links the upload to a registered sequencing VCF
                  (JointCalledVCF / SingleSampleVCF) by path. Defaults to `filename` for backwards
                  compatibility. Pass path=None to omit the query param entirely - required for ad-hoc
                  uploads such as the annotate/download flow, where sending a client-side path makes
                  SeqAuto deployments try (and fail) to match it to a registered VCF.

            metadata: facts about the file it doesn't carry itself, sent as extra query params. A VCF accepts:
                  'genome_build' - the build's own name ('GRCh37'), not an alias ('hg19')
                  'source' - the caller/software, eg 'DRAGEN TSO500 SmallVariant'
                  'extraction' - reference (str or ExternalReference) for every sample in the file
                  'sample_extractions' - {vcf_sample_name: reference} for a multi-sample VCF
                  Send 'extraction' or 'sample_extractions', not both. An unknown key is a 400. An extraction
                  the server doesn't have yet is not - it attaches once the extraction is created.
                  Needs the server feature 'upload_metadata'. Without it, SKIP uploads the file without the
                  metadata (as older clients did) rather than not at all, and ERROR raises

            file_type: the server's name for what this file is, eg UploadFileType.DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT.
                  Not sent (the server decides from the filename) - if given, the upload only happens when
                  accepts_upload(file_type), so an older server doesn't mis-import it as something else """
        if metadata and not self.supports(ServerFeature.UPLOAD_METADATA):
            self._unsupported(f"upload metadata for '{filename}' (server doesn't support feature 'upload_metadata')")
            metadata = None
        if file_type and not self._require_upload(file_type):
            return None
        url = self._get_url("upload/api/v1/file_upload")
        if path is _UNSET:
            path = filename
        params = {"path": path} if path is not None else {}
        if metadata:
            params.update(self._upload_metadata_params(metadata))
        with open(filename, "rb") as f:
            response = requests.post(url, headers=self.headers, files={"file": f}, params=params)
            extra_error_message = f"{filename=}"
            return self._handle_json_response(response, extra_error_message)

    @staticmethod
    def _upload_metadata_params(metadata: dict) -> dict:
        """ Query params are strings - references and objects go as JSON, which the server parses """
        if reserved := {"path", "force"} & set(metadata):
            raise ValueError(f"Upload metadata can't use reserved query param(s): {', '.join(sorted(reserved))}")
        params = {}
        for key, value in metadata.items():
            if isinstance(value, ExternalReference):
                value = reference_json(value)
            elif isinstance(value, dict):
                value = {k: reference_json(v) for k, v in value.items()}
            if isinstance(value, (dict, list)):
                value = json.dumps(value)
            params[key] = value
        return params

    ###############
    ## Get methods

    def _get(self, path, params=None):
        url = self._get_url(path)
        response = requests.get(url,
                                params,
                                headers=self.headers)

        extra_error_message = f"{url=}, {params=}"
        return self._handle_json_response(response, extra_error_message)

    def sequencing_run_has_vcf(self, sequencing_run: SequencingRun, path: Optional[str] = None):
        """ Returns whether the SequencingRun has a VCF associated with it. If path is set, VCF must match that path
            otherwise - any VCF present will retur True """
        return self.sequencing_run_name_has_vcf(sequencing_run.name, path)

    def sequencing_run_name_has_vcf(self, sequencing_run_name: str, path: Optional[str] = None):
        """ Returns whether the SequencingRun has a VCF associated with it. If path is set, VCF must match that path
            otherwise - any VCF present will retur True """
        data = self._get(f"seqauto/api/v1/sequencing_run/{sequencing_run_name}/")
        found_vcf = False
        for vcf in data["vcf_set"]:
            if path is None or path in vcf["path"]:
                found_vcf = True
                break
        return found_vcf

    ##################################
    ## Uploaded file annotation flow

    @staticmethod
    def _upload_key_segment(uploaded_file_id: Optional[int], sha256: Optional[str]) -> str:
        """ URL segment identifying an uploaded file - by uploaded_file_id (preferred) or its SHA-256 hash. """
        if uploaded_file_id is not None:
            return str(uploaded_file_id)
        if sha256:
            return f"sha256/{sha256}"
        raise ValueError("Must provide one of 'uploaded_file_id' or 'sha256'")

    def poll_upload_status(self, uploaded_file_id: Optional[int] = None, sha256: Optional[str] = None) -> Optional[dict]:
        """ Single GET of an uploaded file's import/annotation status.

            Keyed by uploaded_file_id (returned from upload_file) or the SHA-256 of the uploaded file.
            See wait_for_annotation to block until annotation is complete. Needs server feature 'upload_status' """
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        segment = self._upload_key_segment(uploaded_file_id, sha256)
        return self._get(f"upload/api/v1/upload_status/{segment}")

    def wait_for_annotation(self, uploaded_file_id: Optional[int] = None, sha256: Optional[str] = None,
                            timeout: float = 3600, poll_interval: float = 10, sleep: Callable = time.sleep,
                            max_transient_errors: int = 5) -> Optional[dict]:
        """ Poll poll_upload_status until 'annotation_complete' is true, then return the final status dict.

            Raises AnnotationError if the server reports an 'error', or TimeoutError if 'timeout' seconds elapse.
            'sleep' is injectable so tests can avoid real delays.

            Transient server hiccups (5xx / connection / timeout - e.g. a brief 500 right after upload while the
            server is still creating the upload record) are tolerated: up to 'max_transient_errors' *consecutive*
            failures are retried before giving up. A 4xx response is treated as a real error and raised immediately.
            The success counter resets whenever a poll succeeds. Needs server feature 'upload_status' """
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        deadline = time.monotonic() + timeout
        transient_errors = 0
        while True:
            try:
                status = self.poll_upload_status(uploaded_file_id=uploaded_file_id, sha256=sha256)
                transient_errors = 0
            except (requests.HTTPError, requests.ConnectionError, requests.Timeout) as e:
                status_code = getattr(getattr(e, "response", None), "status_code", None)
                if status_code is not None and status_code < 500:
                    raise  # 4xx is a real client error, not a transient blip
                transient_errors += 1
                if transient_errors > max_transient_errors:
                    raise
                self.logger.warning("Transient error polling upload status (%s/%s), retrying: %s",
                                    transient_errors, max_transient_errors, e)
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"Annotation did not complete within {timeout}s (last error: {e})")
                sleep(poll_interval)
                continue
            if error := status.get("error"):
                raise AnnotationError(error, status)
            if status.get("annotation_complete"):
                return status
            if time.monotonic() >= deadline:
                raise TimeoutError(f"Annotation did not complete within {timeout}s (status={status})")
            sleep(poll_interval)

    @staticmethod
    def _attachment_filename(response, export_type: str) -> str:
        content_disposition = response.headers.get("Content-Disposition", "")
        if m := re.search(r'filename\*?=(?:UTF-8\'\')?"?([^";]+)"?', content_disposition):
            return urllib.parse.unquote(m.group(1))
        ext = ".vcf.gz" if export_type == "vcf" else ".csv.zip"
        return f"download{ext}"

    def download_annotated(self, uploaded_file_id: Optional[int] = None, sha256: Optional[str] = None,
                           export_type: str = "vcf", dest_path: Optional[Union[str, Path]] = None,
                           timeout: float = 3600, poll_interval: float = 10,
                           sleep: Callable = time.sleep) -> Optional[Path]:
        """ Download the cohort-level annotated export of an uploaded VCF, saving it to disk.

            export_type is 'vcf' (gzipped *.vcf.gz) or 'csv' (zipped *.csv.zip). The export covers all samples
            in the uploaded VCF (single-sample VCFs included).

            The endpoint returns 202 while the file is still being generated - this polls until it is ready (200),
            then streams the attachment to dest_path. If dest_path is None the server's attachment filename is used
            in the current directory; if it is an existing directory the attachment filename is placed inside it;
            otherwise it is treated as the full destination path. Returns the Path written.

            Raises TimeoutError if the file isn't ready within 'timeout' seconds. Needs server feature 'upload_status' """
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        if export_type not in ("vcf", "csv"):
            raise ValueError(f"export_type must be 'vcf' or 'csv', got {export_type!r}")
        segment = self._upload_key_segment(uploaded_file_id, sha256)
        url = self._get_url(f"upload/api/v1/download/{segment}/{export_type}")
        deadline = time.monotonic() + timeout
        while True:
            response = requests.get(url, headers=self.headers, stream=True)
            if response.status_code == 202:
                if self.log_response:
                    self.logger.info("Download of '%s' still generating: %s", url, self._safe_json(response))
                if time.monotonic() >= deadline:
                    raise TimeoutError(f"Download not ready within {timeout}s (url={url})")
                sleep(poll_interval)
                continue
            if not response.ok:
                # Raises with logging (JSON 'error' message from the server)
                self._handle_json_response(response, f"download {export_type} for {segment}")
            dest = self._resolve_download_dest(dest_path, response, export_type)
            with open(dest, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            return dest

    def _resolve_download_dest(self, dest_path: Optional[Union[str, Path]], response, export_type: str) -> Path:
        filename = self._attachment_filename(response, export_type)
        if dest_path is None:
            return Path(filename)
        dest = Path(dest_path)
        if dest.is_dir():
            return dest / filename
        return dest

    @staticmethod
    def _safe_json(response):
        try:
            return response.json()
        except Exception:
            return None

    def annotate_vcf(self, filename: str, export_type: str = "vcf", dest_path: Optional[Union[str, Path]] = None,
                     timeout: float = 3600, poll_interval: float = 10, sleep: Callable = time.sleep) -> Optional[Path]:
        """ Convenience one-liner: upload a VCF, wait for annotation to finish, download the annotated export.

            Chains upload_file -> wait_for_annotation -> download_annotated and returns the Path written.
            Needs server feature 'upload_status' - checked before uploading """
        if not self._require(ServerFeature.UPLOAD_STATUS):
            return None
        # path is SeqAuto-only and makes ad-hoc uploads fail the import - omit it for the annotate flow
        upload = self.upload_file(filename, path=None)
        uploaded_file_id = upload["uploaded_file_id"]
        self.wait_for_annotation(uploaded_file_id=uploaded_file_id,
                                 timeout=timeout, poll_interval=poll_interval, sleep=sleep)
        return self.download_annotated(uploaded_file_id=uploaded_file_id, export_type=export_type,
                                       dest_path=dest_path, timeout=timeout, poll_interval=poll_interval, sleep=sleep)




