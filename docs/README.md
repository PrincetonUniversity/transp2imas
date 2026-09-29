# transp2imas documentation

This directory contains the maintained technical and developer documentation for
`transp2imas`. The top-level README is intentionally focused on installation and
basic use; the documents here explain the design, validation strategy, mapping
policy, and release process in more detail.

## Start here

- [Architecture](architecture.md) -- module layout, control flow, shared helper
  classes, and where to add new conversion logic.
- [Development guide](development.md) -- local environment, coding conventions,
  Python header policy, and contribution workflow.
- [Testing and validation](testing.md) -- what each `make` target checks and which
  checks require an installed IMAS Data Dictionary.
- [Mapping workflow](mapping-workflow.md) -- the group-by-group procedure used to
  extend TRANSP variable coverage without guessing IMAS semantics.
- [Release process](release-process.md) -- version files, Conventional Commits,
  Release Please, GitHub settings, and how the initial `v0.9.0` baseline is
  established.
- [Technical reference](technical-reference.md) -- detailed CDF-to-IDS mappings,
  units, sign conventions, compatibility notes, provenance, and historical
  development audits.

## Documentation policy

Documentation should be updated in the same pull request as a behavior or
mapping change. In particular:

1. New or changed CDF-to-IMAS mappings belong in `technical-reference.md`.
2. New development or validation procedures belong in `development.md` or
   `testing.md`.
3. Changes to versioning or GitHub Actions belong in `release-process.md`.
4. New IDS groups should be added to the status/roadmap section in
   `mapping-workflow.md`.

The pre-0.9 internal 8.x labels retained in the technical reference are
historical development milestones, not public semantic versions.
