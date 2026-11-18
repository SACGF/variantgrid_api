### Changed

- [`vg_api annotate_vcf --wait` timeout](https://github.com/SACGF/variantgrid_api/issues/52): when `--timeout` runs out before annotation finishes, `vg_api` now says the file is still annotating and to run the same command again later, and exits with 3 (pending). Before, it crashed with a `TimeoutError` traceback. With `--verbose` the last upload status is logged to stderr. No server change needed.
