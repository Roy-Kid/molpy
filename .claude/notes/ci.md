# CI / Pre-commit Conventions

Project-specific parity policy for `/mol:ci-sync` and the `ci-guard` agent.
Migrated from the former local `molpy-ci-sync` skill and `molpy-ci-auditor` agent
(2026-06-10). The authoritative configs are `.pre-commit-config.yaml` and
`.github/workflows/ci.yml` — the pre-commit file header states the parity rule:
every CI check is mirrored as a hook, and any divergence is a bug fixed in the
same commit.

## Canonical check set

| Check | Tool | Pre-commit | CI |
|---|---|---|---|
| Format | `ruff format --check src/ tests/` | required (pre-commit stage) | required |
| Lint | `ruff check src/ tests/` | required (pre-commit stage) | required |
| Type | `ty check src/molpy/` | required (pre-commit stage) | required |
| Tests | `uv run --extra dev python -m pytest tests/ -n auto` | required (**pre-push stage**) | required |
| File hygiene | `pre-commit-hooks` (trailing-ws, eof, merge-conflict, …) | required | intentionally absent |
| Docs build | `zensical build` | intentionally absent (too slow) | required when `docs/` present |

Intentional exemptions: file hygiene hooks are local-only (and skip
`tests/tests-data/`, whose bytes are test inputs);
`zensical build` / docs jobs are CI-only. Do not "fix" these as parity gaps.
Docs deploy is Cloudflare Pages (builds from the repo), not a GitHub workflow.

## Audit rules

- A check counts as "present in CI" only if it runs on `push` to the main branch
  or on `pull_request` — not merely on `workflow_dispatch`.
- Shared tool versions (ruff, ty) must not drift: flag when one side is pinned
  and the other floats, or when major versions differ.
- Ruff settings live in `pyproject.toml` / `ruff.toml`; CI must not override them
  with inline flags (drift risk).
- Install both hook stages: `pre-commit install --hook-type pre-commit --hook-type pre-push`.

## Whole-tree gates

- Lint: `tox -e lint` (ruff format, ruff check, ty check) — the only tox env.
- Tests: `uv run --extra dev python -m pytest tests/ -n auto`, identical in
  ci.yml, full.yml, release.yml and the pre-push hook. It runs in the uv
  project environment because uv honours `[tool.uv.sources]` (the sibling
  molrs checkout); tox's pip does not, which is why the former `tox -e py`
  env was deleted (2026-09-20).

**2026-07-29 incident:** a `language: system` pytest hook ran against an
editable molrs that hid missing PyPI APIs; CI failed while pre-push "passed".
The test hook must always be the CI command above, never a bare `pytest`.

`tox-lint` / `pytest` / `molrs-pin-on-pypi` keep `always_run: true`.

## Ownership split

| Concern | Owner |
|---|---|
| Fixture files | `tests/tests-data/` (committed; read via the `TEST_DATA_DIR` fixture in `tests/conftest.py`) |
| A published `molcrafts-molrs` on the minor line exists on PyPI | `.pre-commit-config.yaml` hook **`molrs-pin-on-pypi`** (pre-push; auto-skips while `[tool.uv.sources]` overrides the pin) |
| Format / lint / type | `tox -e lint` |

## Escape hatch

Hooks can be skipped per-commit without disabling them permanently:

```bash
SKIP=pytest git commit -m "wip: mid-feature"
SKIP=ty,pytest git commit -m "wip: draft"
```

Do **not** `SKIP=molrs-pin-on-pypi` to land features that need unpublished molrs.

## CI matrix convention

**Do not skip a minor** of `requires-python` (`>=3.12` → **3.12, 3.13, 3.14**).

| Workflow | Matrix |
|---|---|
| `ci.yml` `test` | **3 OS** (`ubuntu` / `macos` / `windows`) × **3.12 / 3.13 / 3.14** |
| `full.yml` `test` | same 3×3 on master push |
| `release.yml` | same 3×3 **before** PyPI publish |

Lint stays single-job on ubuntu + 3.14. Never drop 3.13 “because ends only”.
