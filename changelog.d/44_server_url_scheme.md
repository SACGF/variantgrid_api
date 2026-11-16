### Changed

- [`VariantGridAPI(server, ...)` normalises `server`](https://github.com/SACGF/variantgrid_api/issues/44). A server without a scheme (eg `vg_api annotate_vcf --server test.variantgrid.com`, or the same in `VARIANTGRID_API_SERVER`) now defaults to `https://`. Before, it failed with `requests.exceptions.MissingSchema: Invalid URL 'api/v1/capabilities'`. A trailing slash is also added, so a server under a sub-path (`https://example.org/vg`) keeps it instead of `urljoin` dropping it. `VariantGridAPI.server` holds the normalised URL. No server change needed.
