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

There's also a `vg_api annotate_vcf` command line tool (first call uploads, run it again to download once
ready), a step-by-step form (`upload_file` → `wait_for_annotation` → `download_annotated`), and a "submit now,
download later" pattern for long-running jobs. See
**[Annotate a VCF](https://github.com/SACGF/variantgrid_api/wiki/Annotate-a-VCF)** on the wiki.

## Testing

```
# Install required testing packages
python3 -m pip install -e ".[test]"
python3 -m pytest --cov=variantgrid_api
```
