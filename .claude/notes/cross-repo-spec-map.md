# molpy ↔ molrs

**Release order (2026-07-22):** molrs **first**, then molpy. Master landings
for a version need a **tag**. No monorepo merge; no pin-parity scripts.
See `.claude/notes/release.md`.

## Live specs

| molrs | molpy |
|-------|-------|
| (see molrs `.claude/specs/INDEX.md`) | (see molpy `.claude/specs/INDEX.md`) |

## Verbs

molpy's family→verb table lives in molpy `.claude/notes/architecture.md`
§ Design laws 4. Rows and members marked `molrs:` there **describe** molrs and do
not bind it: molpy adopts molrs verbs (`build`, `apply`, `compute`, `typify`,
`generate`, `run`); the reverse never holds (sink direction).

## Version

Minor-line pin in `pyproject.toml`: `molcrafts-molrs>=X.Y.0,<X.(Y+1)` must
include a **published** molrs release on that minor, not only a local rebuild
with a matching version string. Runtime check is major.minor only.
