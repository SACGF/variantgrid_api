## [1.5.0] - Unreleased

### Added

- [Accept sequencing data without FastQs (BAM-first runs)](https://github.com/SACGF/variantgrid_api/issues/18) - `SequencingFile.fastq_r1` / `fastq_r2` are now optional, for sequencers that emit BAM directly or runs where FastQs aren't kept. `create_sequencing_data()` only sends an `unaligned_reads` block when `fastq_r1` is set; otherwise the record is just `bam_file` + `vcf_file` and the server resolves the sample from `sample_name`. Requires the matching server support in SACGF/variantgrid (see SACGF/variantgrid_sapath#357). See `examples/example_bam_first_run.py`.
- [Patient / specimen / extraction, naming a VCF's extraction, and specimen measures](https://github.com/SACGF/variantgrid_api/issues/20) - client support for TSO 500 phase 4. Requires a server at or after SACGF/variantgrid#1716 (SACGF/variantgrid#1707, #1559). See `examples/example_tso500.py`.
  - New dataclasses `Patient`, `Specimen`, `Extraction`, `SpecimenMeasure`, `ExternalPK` and `ExternalReference`, plus enums `Sex`, `TissueStatus`, `NucleicAcid` and `SpecimenMeasureType` holding the values the server stores. Anywhere a parent or extraction is named, pass an `ExternalReference` (local `reference_id` and/or `code` + `external_type` + optional `external_manager`) or a bare string meaning the local reference. `ExternalReference.from_patient()` / `from_specimen()` / `from_extraction()` build one from an object.
  - New client methods `create_patient()`, `create_specimen()`, `create_extraction()`, `create_specimen_measure()`, `create_specimen_measures()` (bulk, one specimen) and `link_sequencing_sample_extraction()`. The creates are upserts, so re-posting a run returns the same rows. The link call returns the server's `match_status` / `match_error`. An extraction the server doesn't have yet comes back as a 202 with `match_status` `Pending` rather than raising, and the link attaches itself once the extraction is created.
  - `upload_file()` takes an optional `metadata` dict, sent as extra query params: `genome_build` (a build's own name such as `GRCh37`, not an alias), `source`, `extraction`, or `sample_extractions` (`{vcf_sample_name: reference}`). References and objects are JSON-encoded. `path` and `force` are reserved and raise `ValueError`.
  - `MockVariantGridAPI` mirrors the new methods and the `metadata` kwarg.
  - `tests/test_data/tso500/` holds the synthetic TSO 500 data set from SACGF/variantgrid.
- [Server capabilities probe and feature-gated calls](https://github.com/SACGF/variantgrid_api/issues/21) - one client can talk to both VG3 and VG4 servers (SACGF/variantgrid_sapath#443). Requires the `seqauto/api/v1/capabilities` endpoint in SACGF/variantgrid; a server without it counts as legacy.
  - `VariantGridAPI.capabilities` is fetched on first use and then cached, and returns a `ServerCapabilities` (`version`, `git_hash`, `features`, `upload_file_types`). A 404 gives `ServerCapabilities.LEGACY`. Building the client and calling ungated methods make no extra request.
  - `supports(feature)` and `accepts_upload(file_type)` let callers branch where skipping isn't enough, such as uploading the CombinedVariantOutput or the splice VCF.
  - New constructor arg `unsupported_feature_policy`: `UnsupportedFeaturePolicy.ERROR` (default) raises `UnsupportedFeatureError`, and `SKIP` logs a warning and returns `None`.
  - Gated calls: `create_patient` / `create_specimen` / `create_extraction` (`patients`), `create_specimen_measure(s)` (`specimen_measures`), `link_sequencing_sample_extraction` (`link_extraction`), and `upload_file` with non-empty `metadata` (`upload_metadata`). For `upload_file`, SKIP drops the metadata but still uploads the file, as older clients did.
  - The upload/annotate flow is now gated on `upload_status`: `poll_upload_status`, `wait_for_annotation`, `download_annotated` and `annotate_vcf` (checked before it uploads). A server without the capabilities endpoint counts as legacy, so these raise `UnsupportedFeatureError` against it even when it has the upload status endpoints. `vg_api annotate_vcf` reports that as an error (exit 1).
  - `upload_file()` takes an optional `file_type` (such as `dragen_tso500_combined_variant_output`). It is not sent to the server; the upload only happens when `accepts_upload(file_type)` is true.
  - `MockVariantGridAPI` takes `capabilities` (default: every gated feature and file type) and `unsupported_feature_policy`. A call skipped under SKIP is not recorded.
  - `examples/example_tso500.py` uses SKIP, and uploads the CombinedVariantOutput when the server accepts it, otherwise the splice VCF.

### Changed

- `SequencingFile` field order is now `sample_name, bam_file, vcf_file, fastq_r1, fastq_r2` (optional fields must come last). Keyword construction is unaffected; positional construction needs updating.
- `upload_file()` without `metadata` sends exactly what it did before.
- `create_sequencing_data()` sends a single-end `unaligned_reads` (no `fastq_r2`) when only `fastq_r1` is set, and raises `ValueError` if `fastq_r2` is set without `fastq_r1`. Records carrying both FastQs send exactly the payload they did before.

## [1.4.0] - 2026-08-07

### Added

- [Joint-called VCFs whose samples span sequencing runs](https://github.com/SACGF/variantgrid_api/issues/19) - `JointCalledVCF` takes an optional `sequencing_samples` list of `SequencingSampleLookup`, so a family trio joint-called from samples sequenced on different runs links all of its members. `sample_sheet_lookup` stays the owning run (the one the VCF path sits under). Requires the matching server support in SACGF/variantgrid, and each contributing run's sample sheet must be sent before the joint call. Asked for in SACGF/variantgrid_sapath#415; see `examples/example_cross_run_trio.py`.

### Changed

- `JointCalledVCF` omits `sequencing_samples` from the posted JSON when it is unset, so a single-run joint call sends exactly the payload it did before and the server keeps deriving members from the owning sample sheet.

## [1.3.2] - 2026-07-13

### Added

- `vg_api annotate_vcf --wait` polls until annotation finishes and downloads it, so you don't need to loop the command yourself.

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
