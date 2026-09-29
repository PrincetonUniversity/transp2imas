# Versioning and release process

## Public versioning

`0.9.0` is the first public pre-1.0 baseline. Earlier 8.x labels in the
technical documentation are internal development snapshots and are not part of
the public semantic-version sequence.

The public version is stored in two places:

- `version.txt` -- canonical Release Please version file;
- `nc2imas.py` -- runtime `__version__`, marked with
  `# x-release-please-version` so Release Please updates it automatically.

`.release-please-manifest.json` records the last established release version.

## Why `.github/` is in the repository

`.github/` is a normal committed repository directory with special meaning to
GitHub. GitHub Actions only discovers workflow YAML files stored under
`.github/workflows/`. Do not put `ci.yml` or `release.yml` at repository root.

This repository uses:

- `.github/workflows/ci.yml` -- portable source checks for pull requests and
  pushes to `main`;
- `.github/workflows/release.yml` -- Release Please automation;
- `.github/PULL_REQUEST_TEMPLATE.md` -- contributor checklist;
- `.github/README.md` -- description of repository-level GitHub automation.

## Initial 0.9.0 bootstrap

Because `0.9.0` is being established from an already-developed code base, treat
it as the existing baseline rather than asking Release Please to infer the old
internal version history.

Recommended bootstrap:

1. merge the 0.9.0 baseline files to `main`;
2. run all available checks;
3. create the Git tag `v0.9.0` and corresponding GitHub Release once;
4. leave `.release-please-manifest.json` at `0.9.0`;
5. from that point forward, let Release Please create subsequent release PRs and
   tags.

## Conventional Commits and version bumps

Release Please derives the next version from conventional commit/PR semantics:

- `fix:` -> patch, for example `0.9.0` to `0.9.1`;
- `feat:` -> minor, for example `0.9.x` to `0.10.0`;
- breaking change -> major.

Documentation/test/chore-only commits normally do not force a user-facing
version bump by themselves.

## Required GitHub repository settings

On GitHub, open:

**Settings -> Actions -> General**

Check the following:

1. GitHub Actions are allowed for the repository under the organization's
   policy.
2. Under **Workflow permissions**, use **Read and write permissions** for the
   release workflow use case.
3. Enable **Allow GitHub Actions to create and approve pull requests** so the
   default `GITHUB_TOKEN` can create the Release Please PR, if organization
   policy permits it.

The source CI workflow itself only needs read access.

## Release lifecycle after 0.9.0

1. Normal feature/fix pull requests are merged to `main`.
2. Release Please updates or opens its release PR.
3. The release PR updates `CHANGELOG.md`, `version.txt`, the manifest, and the
   annotated version in `nc2imas.py`.
4. Review and merge the release PR when ready.
5. Release Please creates the Git tag and GitHub Release.

## If CI should run on Release Please-created PRs

GitHub intentionally limits workflow chaining performed with the repository's
default `GITHUB_TOKEN`. If organization/repository behavior prevents the
Release Please PR from triggering required downstream workflows, use a GitHub
App token or appropriately scoped fine-grained PAT stored as an Actions secret
and pass that token to the Release Please action. Prefer the default token until
this is actually needed.
