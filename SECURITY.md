# Security Policy

## Supported versions

Only the latest released version of AGVarPred is actively supported with
security updates.

| Version | Supported          |
| ------- | ------------------ |
| 1.0.0   | :white_check_mark: |
| < 1.0.0 | :x:                |

## Reporting a vulnerability

If you discover a security issue, please email the maintainers at
<ft.pank@gmail.com> or <kkokay07@gmail.com> rather than opening a public issue.
We will respond as quickly as possible and work with you to assess and address
the issue.

## Known considerations

- The package loads a frozen pickle file (`final_pipeline.pkl`). Always verify
the SHA256 checksum in `manifest.yaml` against a trusted source.
- Never commit API keys or credentials to the repository.
