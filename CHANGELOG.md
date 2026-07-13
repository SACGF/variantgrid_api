## [1.3.1] - 2026-07-13

### Added

- `vg_api` command line tool with an `annotate_vcf` subcommand: first call uploads a VCF, running it again downloads the annotated result once ready.

## [1.3.0] - 2026-07-13

### Added

- [Poll VCF upload status + download annotated VCF/CSV](https://github.com/SACGF/variantgrid_api/issues/17) - New client methods `poll_upload_status()`, `wait_for_annotation()`, `download_annotated()`, and a convenience `annotate_vcf()` wrapper that chains upload → wait → download. Keyed by `uploaded_file_id` or `sha256`. Matches server endpoints added in SACGF/variantgrid#1640. `wait_for_annotation()` tolerates transient 5xx/connection errors while polling (e.g. a brief 500 right after upload).

### Changed

- `upload_file()` now accepts `path=None` to omit the `path` query param. `path` is a SeqAuto-only backend-link hint; sending a client-side path makes SeqAuto deployments try (and fail) to match a registered VCF, so ad-hoc uploads (including `annotate_vcf()`) omit it. Default behaviour is unchanged (still sends `path=filename`).
- [Rename SampleSheetCombinedVCFFile → JointCalledVCF (deprecation)](https://github.com/SACGF/variantgrid_api/issues/16) - New canonical `JointCalledVCF` / `SingleSampleVCF` dataclasses and `create_joint_called_vcf()` client method (POSTs to `seqauto/api/v1/joint_called_vcf/`).

### Deprecated

- `SampleSheetCombinedVCFFile` (use `JointCalledVCF`) and `VCFFile` (use `SingleSampleVCF`) - kept as aliases that emit `DeprecationWarning` on instantiation.
- `create_sample_sheet_combined_vcf_file()` (use `create_joint_called_vcf()`) - kept as a wrapper that emits `DeprecationWarning` and delegates to the new method.

## [1.2.0] - 2026-03-11

### Added

- [Mock VariantGridAPI](https://github.com/SACGF/variantgrid_api/issues/14) - `MockVariantGridAPI` for testing against the client without a live server.
- [Data naming convention helper method](https://github.com/SACGF/variantgrid_api/issues/15)

## [1.1.1] - 2026-01-14

### Added

- [Optional flags to log request / responses](https://github.com/SACGF/variantgrid_api/issues/11)
- [Validation - check for empty data](https://github.com/SACGF/variantgrid_api/issues/12)
- [Add Unit tests](https://github.com/SACGF/variantgrid_api/issues/13)

## [1.1.0] - 2025-10-27

### Added

- [Add datamodels and api calls for Sequencer + SequencerModel](https://github.com/SACGF/variantgrid_api/issues/10)

## [1.0.0] - 2025-08-28

### Changed

- [Single Sample VCF](https://github.com/SACGF/variantgrid_api/issues/7) - Be able to upload upload these as well
- [Make QCExecStats fields optional](https://github.com/SACGF/variantgrid_api/issues/6) - Made optional/mandatory fields match backend

## [0.3.1] - 2025-06-11

### Added

- [New method sequencing_run_has_vcf](https://github.com/SACGF/variantgrid_api/issues/5) - Ability to check whether VCF already associated with sequencing run

### Changed

- [Sequencing Lane should be optional](https://github.com/SACGF/variantgrid_api/issues/2)

## [0.3.0] - 2024-08-22

### Added

- Initial commit in new Repo. Started as v3 as we had very old code on PyPi

[unreleased]: https://github.com/SACGF/variantgrid_api/compare/v1.1.1...HEAD
[1.1.1]: https://github.com/SACGF/variantgrid_api/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/SACGF/variantgrid_api/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/SACGF/variantgrid_api/compare/v0.3.1...v1.0.0
[0.3.1]: https://github.com/SACGF/variantgrid_api/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/SACGF/variantgrid_api/releases/tag/v0.3.0
