# variantgrid_api

[![PyPi version](https://img.shields.io/pypi/v/variantgrid_api.svg)](https://pypi.org/project/variantgrid_api/) [![Python versions](https://img.shields.io/pypi/pyversions/variantgrid_api.svg)](https://pypi.org/project/variantgrid_api/)

Python API client for [VariantGrid](https://github.com/SACGF/variantgrid) Open source Variant database and analysis platform

See [changelog](https://github.com/SACGF/variantgrid_api/blob/main/CHANGELOG.md)

## Install

```
python3 -m pip install variantgrid_api
```

## Example

```
from variantgrid_api.api_client import VariantGridAPI
from variantgrid_api.data_models import EnrichmentKit

api = VariantGridAPI(server="https://variantgrid.com", api_token="YOUR_API_TOKEN")
enrichment_kit = EnrichmentKit(name="idt_haem", version=1)
result = api.create_enrichment_kit(enrichment_kit)
```

## Annotate a VCF and download it back

Upload a VCF, have VariantGrid import + annotate any novel variants, then download the cohort-level annotated
export (all samples, single-sample VCFs included). Annotation can take a while, so the quickest way in is the
`vg_api` command line tool: the first call uploads, and running the same command again downloads the result
once it's ready.

```console
$ export VARIANTGRID_API_TOKEN=YOUR_API_TOKEN
$ vg_api annotate_vcf input.vcf.gz -o results/
Uploaded input.vcf.gz (id=13256).
Annotating input.vcf.gz - run the same command again later to download.

$ vg_api annotate_vcf input.vcf.gz -o results/        # once it's done
Annotated vcf written to results/input.vcf_annotated_v254_GRCh38.vcf.gz
```

From Python it's `upload_file()` / `poll_upload_status()` / `download_annotated()`, or the blocking
`annotate_vcf()` one-liner. See
**[Annotate a VCF](https://github.com/SACGF/variantgrid_api/wiki/Annotate-a-VCF)** on the wiki for batches, the
submit-now/download-later pattern, and all the options.

## Patients, specimens and extractions

VariantGrid can record which patient, specimen and extraction a lab's sequencing came from, and specimen-level
measures such as TMB, MSI and GIS. This needs a server at or after SACGF/variantgrid#1716.

```
from variantgrid_api.data_models import Patient, Specimen, Extraction, ExternalReference, NucleicAcid

api.create_patient(Patient(patient_code="C0000001"))
api.create_specimen(Specimen(patient="C0000001", reference_id="2600000001"))
api.create_extraction(Extraction(specimen="2600000001", reference_id="2600000001C",
                                 nucleic_acid_source=NucleicAcid.DNA))
api.upload_file("sample.vcf.gz", path=None,
                metadata={"extraction": "2600000001C", "genome_build": "GRCh37"})
```

A bare string names a record by its local reference. Use `ExternalReference(code=..., external_type=...)` to
name it by a LIMS identifier instead. See `examples/example_tso500.py` for a full run.

## Talking to more than one VariantGrid version

Servers of different ages accept different calls. The client asks the server which features it has
(`GET seqauto/api/v1/capabilities`, fetched once on first use), and each call that needs a newer server
checks that list first: the patient / specimen / extraction calls, specimen measures,
`link_sequencing_sample_extraction`, the annotate flow (`poll_upload_status`, `wait_for_annotation`,
`download_annotated`, `annotate_vcf`), and `upload_file` with `metadata` or `file_type`. A server without the
capabilities endpoint counts as legacy and reports no features.

By default an unsupported call raises `UnsupportedFeatureError`. To run the same code against old and new
servers, skip those calls instead: they log a warning and return `None`. `upload_file` is the exception: on a
server without upload metadata it drops the `metadata` and still uploads the file.

```
from variantgrid_api.api_client import VariantGridAPI, UnsupportedFeaturePolicy

api = VariantGridAPI(server, api_token, unsupported_feature_policy=UnsupportedFeaturePolicy.SKIP)
api.create_patient(patient)        # skipped on a server without patients

# Uploads only when the server has an importer for this file type
api.upload_file("sample_CombinedVariantOutput.tsv", path=None,
                file_type="dragen_tso500_combined_variant_output")
```

Where the fallback isn't simply "do nothing", branch on `api.supports("feature")` or
`api.accepts_upload("file_type")`. `api.capabilities.version` is useful for logging which server you reached.
For tests, `MockVariantGridAPI(capabilities=ServerCapabilities.LEGACY)` behaves like an old server.

## Testing

```
# Install required testing packages
python3 -m pip install -e ".[test]"
python3 -m pytest --cov=variantgrid_api
```
