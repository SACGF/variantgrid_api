### Changed

- [How to get an API token](https://github.com/SACGF/variantgrid_api/issues/48): VariantGrid users can now create their own token (click your username, then under "API Tokens", "Manage tokens for scripts and API clients"). The README has a new "API token" section, and `vg_api` says how to get a token in `--token`'s help, when no token is given, and when the server rejects one (HTTP 401/403).
