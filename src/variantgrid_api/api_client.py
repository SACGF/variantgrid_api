import datetime
import json
import logging
import urllib
from enum import Enum
from typing import List, Optional, Callable

import requests

from variantgrid_api.data_models import EnrichmentKit, SequencingRun, SampleSheet, SampleSheetCombinedVCFFile, \
    SampleSheetLookup, SequencingFile, QCGeneList, QCExecStats, QCGeneCoverage, SequencerModel, Sequencer


class DateTimeEncoder(json.JSONEncoder):
    def default(self,o):
        if isinstance(o,(datetime.date, datetime.datetime)):
            return o.isoformat()
        return super().default(o)


class EmptyInputPolicy(Enum):
    IGNORE = "ignore"
    WARN = "warn"
    ERROR = "error"


class VariantGridAPI:
    def __init__(self, server, api_token,
                 empty_input_policy=EmptyInputPolicy.ERROR,
                 logger: Optional[logging.Logger] = None,
                 log_request=False, log_response=False
    ):
        self.server = server
        self.headers = {"Authorization": f"Token {api_token}"}
        if logger is None:
            logger = logging.getLogger(__name__)
        self.logger = logger
        self.validation_handler = self._get_validation_handler(empty_input_policy, logger)
        self.log_request = log_request
        self.log_response = log_response

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

    def _validate_list(self, info: str, list_obj: List):
        if not list_obj:
            self.validation_handler(f"{info}: empty list")

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

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file: SampleSheetCombinedVCFFile):
        self._validate_object("sample_sheet_combined_vcf_file", sample_sheet_combined_vcf_file)
        json_data = sample_sheet_combined_vcf_file.to_dict()
        return self._post("seqauto/api/v1/sample_sheet_combined_vcf_file/",
                          json_data)

    def create_sequencing_data(self, sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]):
        self._validate_object("sample_sheet_lookup", sample_sheet_lookup)
        self._validate_list("sequencing_files", sequencing_files)
        records = []
        for sf in sequencing_files:
            data = sf.to_dict()
            # put into hierarchial JSON DRF expects
            fastq_r1 = data.pop("fastq_r1")
            fastq_r2 = data.pop("fastq_r2")
            data["unaligned_reads"] = {
                "fastq_r1": {"path": fastq_r1},
                "fastq_r2": {"path": fastq_r2}
            }
            records.append(data)

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

    def upload_file(self, filename: str):
        url = self._get_url("upload/api/v1/file_upload")
        with open(filename, "rb") as f:
            kwargs = {
                "files": {"file": f},
                "params": {"path": filename}
            }
            response = requests.post(url, headers=self.headers, **kwargs)
            extra_error_message = f"{filename=}"
            return self._handle_json_response(response, extra_error_message)

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




