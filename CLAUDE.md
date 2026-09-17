# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Python client (`variantgrid_api` on PyPI) for the REST API of [VariantGrid](https://github.com/SACGF/variantgrid). The server-side endpoints live in the SACGF/variantgrid repo. Many client features depend on matching server changes, and the CHANGELOG entries link to them.

## Commands

A `.venv` is already set up in the repo.

```bash
python3 -m pip install -e ".[test]"                         # pytest, pytest-cov, responses
python3 -m pytest --cov=variantgrid_api                     # full suite (fast, no network)
python3 -m pytest tests/test_api_client.py::test_create_experiment   # single test
python3 -m build && twine upload dist/*                     # release (dist/ is gitignored)
```

No linter or formatter is configured.

## Architecture

The package is `src/variantgrid_api/` (src layout) and has four modules:

- **`data_models.py`**: `@dataclass_json` dataclasses that mirror the server's SeqAuto models (EnrichmentKit → SequencingRun → SampleSheet → SequencingSample, then JointCalledVCF / SequencingFile (BAM + SingleSampleVCF + optional FastQs) → QC* records). `.to_dict()` produces the JSON the client sends, so the JSON shape is controlled here:
  - Wire names that differ from field names use `field(metadata=config(field_name=...))`. For example, `JointCalledVCF.sample_sheet_lookup` is sent as `sample_sheet`.
  - Optional fields use `config(exclude=lambda x: x is None)`, so leaving a new field unset sends exactly the payload older clients sent. Keep this backwards-compatible pattern when adding fields. Optional fields must come after required ones.
  - Existing records are referenced through `*Lookup` objects (`SampleSheetLookup` = sequencing run name + sample sheet hash, and `SequencingSampleLookup`), not by nesting full objects.
- **`api_client.py`**: `VariantGridAPI`. Each `create_*` method validates its input (`EmptyInputPolicy`: ERROR by default, or WARN/IGNORE), then POSTs to `seqauto/api/v1/...` through `_post`/`_handle_json_response`, which log and `raise_for_status()` on failure. A few methods reshape the output of `to_dict()` into the nested form DRF expects: `create_sample_sheet` replaces `sequencing_run` with its name, and `create_sequencing_data` wraps the FastQs in `unaligned_reads`. The upload/annotation flow (`upload_file` → `poll_upload_status` / `wait_for_annotation` → `download_annotated`, and the `annotate_vcf` wrapper) is under `upload/api/v1/...`. Uploads are keyed by `uploaded_file_id` or `sha256`.
  - **Capability gating** (#21): calls that need a newer server start with `if not self._require("<feature>"): return None` (or `_require_upload(file_type)`), checked against the lazily fetched `capabilities` (`GET api/v1/capabilities`, where a 404 or redirect means `ServerCapabilities.LEGACY`). `unsupported_feature_policy` decides whether that raises `UnsupportedFeatureError` or logs and skips. The exception is `upload_file` metadata: SKIP drops the metadata but still uploads the file. Feature names are the server's contract (`API_FEATURES` in the variantgrid repo's `variantgrid/views_rest.py`). Methods older servers already support stay ungated, so they never trigger the probe.
- **`mock_variantgrid_api.py`**: `MockVariantGridAPI` is a hand-written drop-in that records calls and returns canned values (`set_return`, `get_calls`, `assert_called_once`). **When you add or change a public method on `VariantGridAPI`, update the mock to match** (same signature, deprecation behaviour and capability gating) and add a test in `tests/test_mock_variantgrid_api.py`.
- **`cli.py`**: the `vg_api` console script, which currently has only the `annotate_vcf` subcommand. It keeps no local state, and the server must report the `upload_status` feature. It hashes the VCF and asks the server for status by sha256. On a 404 it uploads; if the file is still pending it exits with 3; once ready it downloads (exit 0). `--wait` polls until done. The token and server come from `--token`/`--server` or `VARIANTGRID_API_TOKEN`/`VARIANTGRID_API_SERVER`.

Gotchas:
- `upload_file(path=...)` is a SeqAuto-only hint that links an upload to a registered VCF. It defaults to the filename for backwards compatibility, but ad-hoc uploads (the annotate flow) must pass `path=None`, or the server's import fails.
- Blocking/polling methods take an injectable `sleep` so tests don't wait in real time.
- Test VCFs sent to a real server need the full contig set in the header, or the server can't detect the genome build.

## Tests

Tests use `responses` to mock HTTP and assert on the posted JSON (`json.loads(responses.calls[-1].request.body)`). `tests/conftest.py` provides the `api` fixture (server `https://example.org`, token `TKN`) and `vg_objects`, a dict containing a complete, realistic object graph for the fake run `Haem_20_999_201231_M02027_0112_000000000_JFT79`, with files under `tests/test_data/idt_haem/`. Extend `vg_objects` rather than building ad-hoc models.

## Deprecations

Renames keep the old name as a deprecated alias that emits `DeprecationWarning`: a subclass with `__post_init__` for dataclasses (`SampleSheetCombinedVCFFile` → `JointCalledVCF`, `VCFFile` → `SingleSampleVCF`), or a wrapper method (`create_sample_sheet_combined_vcf_file`). The mock mirrors these aliases too.

## Examples

`examples/` holds real pipeline integration scripts, not tests. `example_haem_20_999.py` is the full end-to-end ordering of calls; the variant examples are `example_bam_first_run.py` and `example_cross_run_trio.py`. `vg_api.py` and `tau_vg_api.py` are lab pipeline CLIs with extra dependencies (samshee, pandas). New server-facing features usually get an example here. Longer usage docs are on the GitHub wiki, which the README links to.

## Changes and releases

- Every user-visible change goes in `CHANGELOG.md` under the top version (`## [x.y.z] - Unreleased` until released), using Added/Changed/Deprecated sections. Entries link the GitHub issue, name the new methods and fields, and note any server-side dependency.
- Bump `version` in `pyproject.toml` together with the changelog.
- Commit messages reference the issue, e.g. `SACGF/variantgrid_api#18 - <summary>` or `Release 1.3.2: <summary>`.
- `claude/plans/` holds design plans for upcoming work (such as the server capabilities probe).
