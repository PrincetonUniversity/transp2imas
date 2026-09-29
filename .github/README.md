# Repository automation

This directory is intentionally committed to the repository. GitHub uses
`.github/` for repository metadata and automation.

## Contents

- `workflows/ci.yml` runs portable source checks on pull requests and pushes to
  `main`.
- `workflows/release.yml` runs Release Please to maintain the release PR and
  create tags/releases after that PR is merged.
- `PULL_REQUEST_TEMPLATE.md` provides the contributor review checklist.

Workflow files must remain under `.github/workflows/`; GitHub does not discover
Actions workflow YAML stored elsewhere in the repository.

See `docs/release-process.md` for repository settings and the 0.9.0 bootstrap.
