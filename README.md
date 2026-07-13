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

Upload a VCF, wait for VariantGrid to import + annotate any novel variants, then download the
cohort-level annotated export (all samples, single-sample VCFs included). The whole flow is a one-liner:

```python
from variantgrid_api.api_client import VariantGridAPI

api = VariantGridAPI(server="https://variantgrid.com", api_token="YOUR_API_TOKEN")

# export_type is "vcf" (gzipped *.vcf.gz) or "csv" (zipped *.csv.zip)
path = api.annotate_vcf("input.vcf", export_type="vcf", dest_path="/data/results/")
print(f"Annotated VCF written to {path}")
```

Or drive each step yourself:

```python
# For ad-hoc uploads pass path=None so the upload isn't treated as a SeqAuto backend-link hint
upload = api.upload_file("input.vcf", path=None)
uploaded_file_id = upload["uploaded_file_id"]           # or upload["sha256_hash"]
api.wait_for_annotation(uploaded_file_id, timeout=3600, poll_interval=10)  # raises on error/timeout
path = api.download_annotated(uploaded_file_id, export_type="csv", dest_path="/data/results/")
```

`poll_upload_status(uploaded_file_id)` returns the raw status dict (including `annotation_complete`,
`progress_percent`, `error`, `vcf_id`, `samples`, ...) if you want to inspect progress directly. All of
these accept `sha256=<hash>` instead of `uploaded_file_id` - the server dedups on the content hash, so it is
stable across machines and useful if you didn't retain the id.

Notes:

- The `annotate_vcf` wrapper already uploads with `path=None`. Only pass a `path` to `upload_file` for SeqAuto
  uploads that link to a registered `JointCalledVCF` / `SingleSampleVCF` (that is what `path` is for - it is
  ignored on non-SeqAuto deployments).
- Downloads require the target deployment to have the cohort export analysis templates configured
  (`ANALYSIS_TEMPLATES_AUTO_COHORT_EXPORT`) - otherwise a clear error is returned.

## Testing

```
# Install required testing packages
python3 -m pip install -e ".[test]"
python3 -m pytest --cov=variantgrid_api
```
