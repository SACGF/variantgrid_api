import datetime
import json
import logging
import os
from typing import List, Optional

import requests

from variantgrid_api.data_models import EnrichmentKit, SequencingRun, SampleSheet, SampleSheetCombinedVCFFile, \
    SampleSheetLookup, SequencingFile, QCGeneList, QCExecStats, QCGeneCoverage, SequencerModel, Sequencer


class DateTimeEncoder(json.JSONEncoder):
    def default(self,o):
        if isinstance(o,(datetime.date, datetime.datetime)):
            return o.isoformat()
        return super().default(o)



class VariantGridAPI:
    def __init__(self, server, api_token):
        self.server = server
        self.headers = {"Authorization": f"Token {api_token}"}

    def _get_url(self, url):
        return os.path.join(self.server, url)

    def _post(self, path, json_data):
        url = self._get_url(path)
        json_string = json.dumps(json_data, cls=DateTimeEncoder)
        response = requests.post(url,
                                 headers={**self.headers, "Content-Type": "application/json"},
                                 data=json_string)

        extra_error_message = f"{url=}, {json_string=}"
        return self._handle_json_response(response, extra_error_message)

    def _handle_json_response(self, response, extra_error_message: Optional[str] = None):
        try:
            json_response = response.json()
        except Exception as e:
            json_response = f"Couldn't convert JSON: {e}"
        if not response.ok:
            if extra_error_message:
                logging.error(extra_error_message)
            logging.error("Response: %s", json_response)
            response.raise_for_status()
        return json_response

    def create_experiment(self, experiment: str):
        json_data = {
            "name": experiment
        }
        return self._post("seqauto/api/v1/experiment/", json_data)

    def create_enrichment_kit(self, enrichment_kit: EnrichmentKit):
        return self._post("seqauto/api/v1/enrichment_kit/",
                          enrichment_kit.to_dict())


    def create_sequencer_model(self, sequencer_model: SequencerModel):
        return self._post("seqauto/api/v1/sequencer_model/",
                          sequencer_model.to_dict())

    def create_sequencer(self, sequencer: Sequencer):
        return self._post("seqauto/api/v1/sequencer/",
                          sequencer.to_dict())

    def create_sequencing_run(self, sequencing_run: SequencingRun):
        return self._post("seqauto/api/v1/sequencing_run/",
                          sequencing_run.to_dict())

    def create_sample_sheet(self, sample_sheet: SampleSheet):
        json_data = sample_sheet.to_dict()
        # We don't want all sequencing_run just the name
        sequencing_run = json_data.pop("sequencing_run")
        json_data["sequencing_run"] = sequencing_run["name"]
        return self._post("seqauto/api/v1/sample_sheet/",
                          json_data)

    def create_sample_sheet_combined_vcf_file(self, sample_sheet_combined_vcf_file: SampleSheetCombinedVCFFile):
        json_data = sample_sheet_combined_vcf_file.to_dict()
        return self._post("seqauto/api/v1/sample_sheet_combined_vcf_file/",
                          json_data)

    def create_sequencing_data(self, sample_sheet_lookup: SampleSheetLookup, sequencing_files: List[SequencingFile]):
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
        json_data = qc_gene_list.to_dict()
        return self._post("seqauto/api/v1/qc_gene_list/",
                          json_data)


    def create_multiple_qc_gene_lists(self, qc_gene_lists: List[QCGeneList]):
        json_data = {
            "records": [
                qcgl.to_dict() for qcgl in qc_gene_lists
            ]
        }
        return self._post("seqauto/api/v1/qc_gene_list/bulk_create",
                          json_data)

    def create_qc_exec_stats(self, qc_exec_stats: QCExecStats):
        json_data = qc_exec_stats.to_dict()
        return self._post("seqauto/api/v1/qc_exec_summary/",
                          json_data)

    def create_multiple_qc_exec_stats(self, qc_exec_stats: List[QCExecStats]):
        json_data = {
            "records": [
                qces.to_dict() for qces in qc_exec_stats
            ]
        }
        return self._post("seqauto/api/v1/qc_exec_summary/bulk_create",
                          json_data)

    def create_multiple_qc_gene_coverage(self, qc_gene_coverage_list: List[QCGeneCoverage]):
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




