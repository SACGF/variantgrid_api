### Added

- [`vg_api annotate_vcf --verbose`](https://github.com/SACGF/variantgrid_api/issues/25) (`-v`) - logs each request and response, and the full upload status JSON on the pending and error paths, to stderr.

### Changed

- [`vg_api annotate_vcf` says why a file isn't ready](https://github.com/SACGF/variantgrid_api/issues/25) - the pending message names the gate that hasn't passed, from the upload status (`pipeline_status`, `import_status` / `progress_percent`, `remaining_annotation_runs`, `annotation_complete`, `downloads_available`), eg `not ready yet - pipeline Success, import 100.0%, 2 annotation runs remaining`. Before, it showed only the import progress, which reads 100% while annotation is still running. Exit codes are unchanged.
