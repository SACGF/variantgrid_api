### Added

- [`ServerFeature.SYMBOLIC_ALTS`](https://github.com/SACGF/variantgrid_api/issues/38) (`"symbolic_alts"`) - the server imports VCFs with symbolic alts (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`) as structural variants. VG4 reports it from SACGF/variantgrid#2732. VG3 never does.
  - `upload_file(..., required_features=[...])` - server features the file needs, checked before uploading like `file_type`. Without them `unsupported_feature_policy` applies: ERROR raises `UnsupportedFeatureError`, SKIP logs and doesn't upload. Nothing extra is sent. Use `required_features=[ServerFeature.SYMBOLIC_ALTS]` for an SV/CNV VCF with symbolic alts. Gene-level files keep their own `file_type` check (`UploadFileType.GENE_LEVEL_CNV_VCF`, `UploadFileType.GENE_LEVEL_SPLICE_VCF`).
  - `MockVariantGridAPI.upload_file()` takes `required_features`, gates on it the same way, and records it only when given.
  - `examples/example_tso500.py` uploads the `_DragenExonCNV.vcf` (coordinate `<DEL>` calls) only to a server with `symbolic_alts`, and the SEGID `cnv.vcf` only to a server that accepts `gene_level_cnv_vcf`.

### Changed

- [`VariantGridAPI(server, api_token=None)` sends no `Authorization` header](https://github.com/SACGF/variantgrid_api/issues/38) - enough to read `capabilities` from a server at or after SACGF/variantgrid#2732, which makes the endpoint public. Before, the client sent `Token None`, which the server rejects. Anonymous callers get no `git_hash` (already optional on `ServerCapabilities`) and a `version` without the commit suffix. Every other call still needs a token.
- README: removed the stale mention of specimen measures (removed in 1.9.0) from the capability-gated calls.
