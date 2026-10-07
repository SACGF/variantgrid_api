### Added

- [`UploadFileType.GENE_LEVEL_SPLICE_VCF`](https://github.com/SACGF/variantgrid_api/issues/30) (`"gene_level_splice_vcf"`) - SpliceGirl's TSO 500 `SpliceVariants.vcf`, imported as gene-level splice variants by a server at or after SACGF/variantgrid#1903. `examples/example_tso500.py` uses it to tell whether the server takes the splice calls from that VCF or (older servers) from the CombinedVariantOutput.

### Changed

- [DRAGEN TSO 500 run-level uploads take `sequencing_run` metadata](https://github.com/SACGF/variantgrid_api/issues/30) - documented in the README, in `upload_file()` and on `UploadFileType`. `DRAGEN_TSO500_METRICS_OUTPUT` requires it. For `DRAGEN_TSO500_COMBINED_VARIANT_OUTPUT` the file names its run 'NA', so the server takes the run from this metadata; without it, it uses the registered run whose current sample sheet names the pair's sample IDs, and the import fails if there's none (SACGF/variantgrid#1904). It's the only metadata either file takes - a current server rejects `genome_build` on a CombinedVariantOutput.
- `examples/example_tso500.py` no longer parses TMB / MSI / GIS out of the CombinedVariantOutput and posts them as specimen measures. It uploads the CombinedVariantOutput and MetricsOutput with `sequencing_run` metadata, and the splice VCF whenever the server has a splice VCF importer.
- README: SeqAuto writes (the sequencing `create_*` calls, the QC calls and `link_sequencing_sample_extraction`) need a superuser, or a member of the server's SeqAuto write group (`seqauto_api_write` by default), otherwise the server returns 403. Reads are unchanged.

### Deprecated

- [Specimen measures](https://github.com/SACGF/variantgrid_api/issues/30) - the server removed them (SACGF/variantgrid#1904): TMB / MSI / GIS now come from the uploaded CombinedVariantOutput, and a current server no longer reports `specimen_measures`. Each of these warns (`DeprecationWarning`) and still works against an older server that reports the feature; against a current one the calls stay gated as before (`UnsupportedFeatureError`, or a skip under `UnsupportedFeaturePolicy.SKIP`), never reaching the removed URLs:
  - `VariantGridAPI.create_specimen_measure()` and `create_specimen_measures()`, and the same on `MockVariantGridAPI` (whose default capabilities still include `specimen_measures`, so the calls still record).
  - `SpecimenMeasure` warns when constructed, `SpecimenMeasureType` when a member is named (`SpecimenMeasureType.TMB`). JSON is unchanged.
  - `ServerFeature.SPECIMEN_MEASURES` - kept, but `supports()` (real and mock) warns when asked for it. Deprecated features are listed in `data_models.DEPRECATED_SERVER_FEATURES`.
