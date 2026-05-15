# h5cpp AI Skill Knowledge Base

Curated, test-verified patterns for AI coding assistants. Every snippet tagged
`tested: true` appears verbatim in a ctest. Run `scripts/verify-skills` in CI
to detect drift between skill files and the test suite.

## How to use

**Claude / Kimi / Codex:**
```bash
scripts/mem search "h5cpp <topic>"   # before starting any h5cpp task
```

**Copilot:** reads this directory via `.github/copilot-instructions.md`.

## Skill files

| File | Topics | Tags |
|---|---|---|
| [includes.md](includes.md) | entry points, include order, h5cpp/all vs core+io | `includes` |
| [io-create.md](io-create.md) | file and dataset creation, RAII handles | `io` `create` `file` `dataset` |
| [io-write.md](io-write.md) | full write, partial write, offset/stride | `io` `write` `partial` |
| [io-read.md](io-read.md) | full read, partial read, type constraints | `io` `read` `partial` |
| [io-append.md](io-append.md) | packet table, streaming append, flush | `io` `append` `streaming` `pt` |
| [properties.md](properties.md) | property list composition, `operator\|`, dcpl/fapl/lcpl | `properties` `dcpl` `chunk` `compression` |

## Tag index

| Tag | Skills |
|---|---|
| `includes` | includes.md |
| `io` | io-create, io-write, io-read, io-append |
| `create` | io-create |
| `write` | io-write |
| `read` | io-read |
| `append` `streaming` `pt` | io-append |
| `chunk` `compression` `gzip` | properties, io-create |
| `partial` `offset` `stride` | io-write, io-read |
| `stl` `vector` | io-write, io-read |
| `string` | io-write, io-read |
| `properties` `dcpl` `fapl` `lcpl` | properties |

## Constraints and known limitations

- `h5::read<std::string>` does **not** compile — use `h5::read<std::vector<std::string>>`.
- True scalar dataset round-trip (non-vector) is blocked on issue #89. Use `std::vector<T>{value}`.
- All snippets require C++17.
