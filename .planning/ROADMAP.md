# Roadmap: HCVTyper — De novo-Informed Strain Selection

## Milestones

- ✅ **v1.0 De novo Minor Confirmation** — Phases 1-4 (shipped 2026-06-08) — see [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)
- ⏸️ **v2.0 HCVGLUE Refactor** — Phase 5 (paused mid-execution, 2/3 plans — see `.continue-here.md`)
- ✅ **v3.0 Strain Model Redraw** — Phases 6-9 (shipped 2026-06-14) — see [milestones/v3.0-ROADMAP.md](milestones/v3.0-ROADMAP.md)

## Phases

<details>
<summary>✅ v1.0 De novo Minor Confirmation (Phases 1-4) — SHIPPED 2026-06-08</summary>

- [x] Phase 1: Major-Gate + Selection-Script Hardening (3/3 plans) — completed 2026-06-06
- [x] Phase 2: De Novo Evidence Plumbing (3/3 plans) — completed 2026-06-07
- [x] Phase 3: De Novo Confirmation Logic (3/3 plans) — completed 2026-06-07
- [x] Phase 4: Regression Guard + Reproducibility (3/3 plans) — completed 2026-06-07

Full detail: [milestones/v1.0-ROADMAP.md](milestones/v1.0-ROADMAP.md)

</details>

### v2.0 HCVGLUE Refactor (PAUSED)

- [ ] **Phase 5: HCVGLUE Per-Sample Parallel Refactor** — PAUSED at Plan 03 Task 2 (blocked by `run-glue.sh` container defects, see `.continue-here.md`). Resume via `/gsd-resume-work`.

<details>
<summary>✅ v3.0 Strain Model Redraw (Phases 6-9) — SHIPPED 2026-06-14</summary>

- [x] Phase 6: Neutral Candidate Selection (3/3 plans) — completed 2026-06-12
- [x] Phase 7: Per-Genotype Assembly Support (2/2 plans) — completed 2026-06-13
- [x] Phase 8: Dominance Scoring + Strain-Role Classification (2/2 plans) — completed 2026-06-14
- [x] Phase 9: Compatibility, Filename Migration + Regression Suite (4/4 plans) — completed 2026-06-14

Full detail: [milestones/v3.0-ROADMAP.md](milestones/v3.0-ROADMAP.md)

</details>

## Progress

| Phase | Milestone | Plans Complete | Status | Completed |
| ----- | --------- | -------------- | ------ | --------- |
| 1. Major-Gate + Selection-Script Hardening | v1.0 | 3/3 | Complete | 2026-06-06 |
| 2. De Novo Evidence Plumbing | v1.0 | 3/3 | Complete | 2026-06-07 |
| 3. De Novo Confirmation Logic | v1.0 | 3/3 | Complete | 2026-06-07 |
| 4. Regression Guard + Reproducibility | v1.0 | 3/3 | Complete | 2026-06-07 |
| 5. HCVGLUE Per-Sample Parallel Refactor | v2.0 | 2/3 | Paused | — |
| 6. Neutral Candidate Selection | v3.0 | 3/3 | Complete | 2026-06-12 |
| 7. Per-Genotype Assembly Support | v3.0 | 2/2 | Complete | 2026-06-13 |
| 8. Dominance Scoring + Strain-Role Classification | v3.0 | 2/2 | Complete | 2026-06-14 |
| 9. Compatibility, Filename Migration + Regression Suite | v3.0 | 4/4 | Complete | 2026-06-14 |

## Next Milestone

Post-v3.0 candidates (see PROJECT.md):

- Finish v2.0 HCVGLUE refactor (Phase 5, paused) — resume via `/gsd-resume-work`
- Close Nyquist validation (phases 01–03) + fix the latent Phase-1 `no_mapping` crash
- Drop the legacy `Major_*`/`Minor_*` column aliases after the one-release deprecation window (COMPAT-03)
- Run Docker test suite confirmation for Phase 8 (tidyverse container)
- Regenerate workflow nf-test snapshot (ENOSPC deferred from v3.0)
