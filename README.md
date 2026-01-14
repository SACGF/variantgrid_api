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

## Testing

```
# Install required testing packages
python3 -m pip install -e ".[test]"
python3 -m pytest --cov=variantgrid_api
```
