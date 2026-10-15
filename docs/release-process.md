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

## `feat:` troubleshooting and pre-1.0 release semantics

The repository explicitly configures the `default` versioning strategy and
`bump-patch-for-minor-pre-major: false`. Accordingly a properly parsed
`feat:` commit requests a **minor** bump even on `0.x`, while `fix:` requests
a patch bump. Both prefixes have visible changelog sections; `perf`, `docs`,
`test`, `refactor` and `chore` are classified separately. Only qualifying
releasable changes cause a version bump; `docs:` and `chore:` alone normally
do not. These settings are verified by `python examples/verify_release_config.py`.

Release Please acts on **commits actually merged into `main`**, not on a PR's
internal commits if they were replaced by a squash merge subject. Set the PR
**title/squash commit subject** to the intended Conventional Commit syntax:

```text
feat(theory): support GLF23 and MMM07 in core_transport
fix(imas): correct TIME3 dimensions
perf: reduce CDF rereads
```

Do not use `Feature: ...`, `Feat: ...`, or only a generic
`Merge pull request #...` subject; they are not Conventional Commit headers.
If `fix:` creates a release PR but `feat:` does not, inspect **the actual
merged commit subject**, Release Please action log, current release PR and
manifest/tag alignment, instead of assuming that adding `feat` to a changelog
section corrects GitHub history:

```bash
git log --first-parent -n 15 --format='%h %s' main
git tag --list 'v*' --sort=-v:refname | head
cat .release-please-manifest.json
```

`bootstrap-sha` applies to the **first** manifest release window; if the
repository has since released through Release Please, the last merged release
PR determines the subsequent commit range. Do not reset `bootstrap-sha`,
manifest or tags to force unrecognized feature commits without reviewing
that history. No local file-only test can prove how your remote release job
handled a past merge; inspect the live workflow log if symptoms persist.

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
