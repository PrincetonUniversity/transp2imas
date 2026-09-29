## Summary

Describe what changes and why.

## Mapping/physics impact

- [ ] No physics/data mapping changed.
- [ ] Mapping changed; source quantity, units, grid, sign convention, and IMAS
      destination are documented.
- [ ] A source-side validation or regression check was added/updated where
      practical.

## Compatibility

- [ ] DD3 behavior considered where supported.
- [ ] DD4/COCOS17 behavior considered where supported.
- [ ] Missing optional source data is handled without fabricating values.

## Validation

- [ ] `make check`
- [ ] `make verify-all`
- [ ] `make verify-dd` (when applicable and an IMAS environment is available)

## Documentation

- [ ] README/docs updated for user-visible or developer-visible behavior.
- [ ] Technical mapping reference updated for mapping changes.
- [ ] Pull-request title follows Conventional Commits (`fix:`, `feat:`,
      `docs:`, `test:`, etc.).
