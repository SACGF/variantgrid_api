### Changed

- [`download_annotated()` / `annotate_vcf()` with `compression="none"` under `UnsupportedFeaturePolicy.ERROR`](https://github.com/SACGF/variantgrid_api/issues/43): the `UnsupportedFeatureError` message no longer ends with "- downloading the zip", because nothing is downloaded. The SKIP warning still says it. `MockVariantGridAPI` matches.
- [`download_annotated(dest_path=Path(...))`](https://github.com/SACGF/variantgrid_api/issues/43) logs a warning when the `Path` doesn't exist and has no suffix (eg `Path("out/new_dir/")`, which `Path` stores as `out/new_dir`). It is still written as a file, as before. To download into a new directory, pass a str ending in `/`.
- Removed the leftover `changelog.d/38_symbolic_alts.md` fragment, already merged into `[1.11.0]`.
