# Branch and PR Governance

## Branch naming

All work branches follow the pattern:

```
<issue-number>-<slug>
```

Examples: `82-modernize-cmake`, `87-datatype-synthesis`, `91-fix-append-bugs`

The legacy `issue/<number>-<slug>` prefix is **rejected** by `new-worktree.sh` and must not be used.

## Creating a worktree

Every issue gets its own git worktree. Never work directly in the canonical repo at
`/home/steven/projects/h5cpp`. Use the workspace script:

```bash
# From /home/steven/projects/vargalabs-workspace

# Preferred: let gh derive the branch name from the issue title
bash scripts/new-worktree.sh h5cpp --issue 99

# With an explicit slug
bash scripts/new-worktree.sh h5cpp --issue 99 --slug my-slug

# Explicit branch name (must match ^[0-9]+-[a-z0-9][a-z0-9-]*$)
bash scripts/new-worktree.sh h5cpp 99-my-slug
```

Worktrees land in:

```
/home/steven/projects/vargalabs-workspace/worktrees/h5cpp/<branch>/
```

The default base is `origin/release`. Pass an explicit ref as the last argument to
override (e.g. `origin/87-datatype-synthesis` for a donor-dependent branch).

## Branch types

### Regular issue branch

Based on `origin/release`. This is the default. PR targets `release`.

### Donor-dependent branch

Based on another issue branch rather than `release`. Use when the issue's work
depends on changes that are not yet merged.

```bash
bash scripts/new-worktree.sh h5cpp --issue 85 --slug detection-idiom \
    origin/87-datatype-synthesis
```

**Rules for donor-dependent branches:**

1. Do not open a PR until the donor branch is merged to `release`.
2. After the donor merges, rebase the dependent on the updated `release` before
   opening the PR.
3. Document the dependency in the PR description (`Depends on #N`).

Current example: `85-refactor-normalize-detection-idiom-usage-across-meta-layer`
is based on `87-datatype-synthesis`.

## Merge targets and flow

```
<number>-<slug>  →  staging  →  release
```

Issue branches PR into `staging`. Once the staging CI matrix is green and the
batch is reviewed, `staging` is merged into `release` for the official cut.

| Branch | Purpose | PR target |
|--------|---------|-----------|
| `release` | Stable, tagged releases | — |
| `staging` | Integration / pre-release validation | `release` |
| `<number>-<slug>` | Issue work | `staging` |

CI runs on push to `release` and `staging`, and on all PRs targeting those branches.

## PR conventions

- **Title:** `[#N] short description` — matches the commit prefix style
- **Body:** include `Closes #N` so GitHub links and auto-closes the issue on merge
- **One issue per PR.** Split unrelated fixes into separate branches.

## Commit message prefix

```
[#N]:author: verb, description
```

Example: `[#91]:amelia: fix aligned_alloc UB in pipeline_t`

## Parallel work

Independent issues (no donor relationship) can be worked simultaneously on separate
worktrees. The workspace supports up to one worktree per branch. Issues that share a
donor must be serialised: donor first, then dependents in order.

## Closing issues

- Issues resolved by a PR are closed automatically via `Closes #N` on merge.
- Issues fixed as a side effect of another PR (e.g. community bug reports resolved
  by a refactor) should be closed manually with a comment citing the PR.
- Issues determined to be already complete should be closed with a comment
  explaining what delivered the fix.

## Worktree cleanup

After a PR merges, remove the local worktree:

```bash
cd /home/steven/projects/h5cpp
git worktree remove \
    /home/steven/projects/vargalabs-workspace/worktrees/h5cpp/<branch>
git branch -d <branch>
```
