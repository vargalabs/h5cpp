# H5CPP Workflow

This document defines the issue-based development workflow for H5CPP.

## Branch Naming Convention

All active development work must use the following branch naming pattern:

```
issue/<number>-<slug>
```

Examples:
- `issue/82-modernize-cmake`
- `issue/95-doctest`
- `issue/86-type-families`

## Donor Branch Policy

Donor branches are read-only archaeological references. They are mined for commits,
never merged into `master` directly.

| Donor Branch | Feeds Issues |
|--------------|-------------|
| `82-build-review-and-modernize-cmake-for-c++17-and-package-exports` | #82, #94, #95 |
| `84-build-add-github-pages-documentation-site-and-ci-integration` | #84 |
| `wip/laptop-typesystem-import-20260507` | #85–#90 |
| `rework` | #85–#90 |
| `dev`, `staging` | Future API features (post-#93) |
| `debian` | Future packaging issue |
| `threaded_compressors` | Future feature |

## Archive Tags

The following branches have been retired from active validation and preserved as
read-only archive tags:

- `archive/appveyor` (was `AppVeyor`)
- `archive/ubuntu-1804` (was `ubuntu-18.04LTS`)
- `archive/vs2017` (was `vs2017-windows`)

## Master Protection Rules

- **Require linear history** — rebase merge only; no merge bubbles
- **Require pull request review** — 1 approval minimum before merging
- **No direct pushes to `master`** — all changes via PR
- **Require status checks** — enforced once CI lands via #84

## Execution Order

Canonical modernization sequence:

```
#82 → #95 → #94 → #85 → #86 → #88 → #90 → #87 → #89 → #91/#92/#93
```

**Hard rule:** No `h5cpp/H5Tmeta.hpp` work before #95 and #94 Phase A are merged.
