# Screen results — candidates the model cannot rule out

Generated 2026-09-26 by `PYTHONPATH=. python -m screen.report`. 98 candidates: every combination of up to 2 of the model's 7 named intervention points, each turned both ways at potency 0.5.

**This is a kill list with a remainder, not a ranking.** Survivors are listed ALPHABETICALLY, never by predicted benefit. Ranking is gated on backtest/lomo.py beating its null; it does not (45.4pp vs 11.8pp).

## Counts

| verdict | n |
|---|---|
| DEGENERATE | 12 |
| OUT_OF_REGIME | 32 |
| SURVIVED | 37 |
| UNREACHABLE | 17 |

## Survivors (37)

Alphabetical. The damage column is a simulated median over 8 infection histories, shown so a reader can challenge the verdict — not as a score.

| candidate | best simulated damage | vs untreated | at potency |
|---|---|---|---|
| `alpha_E+,gamma_R-` | 0.0102 | -99% | 0.8 |
| `alpha_E-,delta+` | 0.2951 | -80% | 0.8 |
| `alpha_E-,delta-` | 0.0573 | -96% | 0.8 |
| `alpha_E-,gamma_E+` | 0.2917 | -80% | 0.8 |
| `alpha_E-,gamma_E-` | 0.0694 | -95% | 0.8 |
| `alpha_E-,gamma_R-` | 0.0014 | -100% | 0.8 |
| `alpha_E-,ke-` | 0.0076 | -99% | 0.8 |
| `alpha_E-,naive_E+` | 0.2333 | -84% | 0.8 |
| `alpha_E-,naive_E-` | 0.5195 | -64% | 0.8 |
| `alpha_R+,gamma_R+` | 1.2214 | -16% | 0.8 |
| `alpha_R+,gamma_R-` | 0.0017 | -100% | 0.8 |
| `alpha_R+,ke+` | 0.1636 | -89% | 0.8 |
| `alpha_R+,ke-` | 0.0028 | -100% | 0.8 |
| `alpha_R+,naive_E+` | 0.0407 | -97% | 0.8 |
| `alpha_R+,naive_E-` | 0.0591 | -96% | 0.8 |
| `delta+,ke-` | 0.3056 | -79% | 0.8 |
| `delta-,ke-` | 0.1037 | -93% | 0.8 |
| `gamma_E+,alpha_R+` | 0.0616 | -96% | 0.8 |
| `gamma_E+,delta-` | 0.1542 | -89% | 0.8 |
| `gamma_E+,gamma_R-` | 0.0011 | -100% | 0.8 |
| `gamma_E+,ke-` | 0.3677 | -75% | 0.8 |
| `gamma_E-,alpha_R+` | 0.0491 | -97% | 0.8 |
| `gamma_E-,delta+` | 0.4901 | -66% | as enumerated |
| `gamma_E-,delta-` | 0.6548 | -55% | as enumerated |
| `gamma_E-,gamma_R-` | 0.0314 | -98% | 0.8 |
| `gamma_E-,ke-` | 0.2473 | -83% | as enumerated |
| `gamma_E-,naive_E+` | 0.5692 | -61% | as enumerated |
| `gamma_E-,naive_E-` | 0.6388 | -56% | 0.2 |
| `gamma_R-` | 0.0027 | -100% | 0.8 |
| `gamma_R-,delta+` | 0.0027 | -100% | 0.8 |
| `gamma_R-,delta-` | 0.0026 | -100% | 0.8 |
| `gamma_R-,ke+` | 0.0053 | -100% | 0.8 |
| `gamma_R-,ke-` | 0.0021 | -100% | 0.8 |
| `gamma_R-,naive_E+` | 0.0052 | -100% | 0.8 |
| `gamma_R-,naive_E-` | 0.0014 | -100% | 0.8 |
| `naive_E+,ke-` | 0.0597 | -96% | 0.8 |
| `naive_E-,ke-` | 0.2812 | -81% | 0.8 |

## Killed

| candidate | verdict | why |
|---|---|---|
| `alpha_E+,alpha_R+` | DEGENERATE | same intervention points as lenercept |
| `alpha_E-` | DEGENERATE | same intervention points as ustekinumab |
| `alpha_E-,alpha_R+` | DEGENERATE | same intervention points as lenercept |
| `alpha_R+` | DEGENERATE | same intervention points as daclizumab |
| `alpha_R+,delta+` | DEGENERATE | same intervention points as glatiramer acetate |
| `alpha_R+,delta-` | DEGENERATE | same intervention points as glatiramer acetate |
| `delta+,naive_E+` | DEGENERATE | same intervention points as IFN-gamma |
| `delta-` | DEGENERATE | same intervention points as abatacept |
| `delta-,naive_E+` | DEGENERATE | same intervention points as IFN-gamma |
| `delta-,naive_E-` | DEGENERATE | same intervention points as IFN-gamma |
| `gamma_E-` | DEGENERATE | same intervention points as Tovaxin |
| `ke-` | DEGENERATE | same intervention points as rituximab |
| `alpha_E+,alpha_R-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E+,gamma_E+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E+,gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E+,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E-,alpha_R-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E-,gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,delta+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,delta-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,ke-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,naive_E+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_R-,naive_E-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `delta+,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,alpha_R-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,delta+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,naive_E+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E+,naive_E-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E-,alpha_R-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_E-,gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,delta+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,delta-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,ke-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,naive_E+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `gamma_R+,naive_E-` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `naive_E-,ke+` | OUT_OF_REGIME | at least one infection history diverged; damage is undefined |
| `alpha_E+` | UNREACHABLE | best probed damage 4.4073 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,delta+` | UNREACHABLE | best probed damage 8.7195 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,delta-` | UNREACHABLE | best probed damage 6.9485 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,gamma_E-` | UNREACHABLE | best probed damage 5.9400 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,ke-` | UNREACHABLE | best probed damage 1.6474 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,naive_E+` | UNREACHABLE | best probed damage 4.7360 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,naive_E-` | UNREACHABLE | best probed damage 51.4810 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,ke+` | UNREACHABLE | best probed damage 1.4137 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,gamma_R-` | UNREACHABLE | best probed damage 3.5176 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `delta+` | UNREACHABLE | best probed damage 3.4541 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `delta+,naive_E-` | UNREACHABLE | best probed damage 11.7083 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `delta-,ke+` | UNREACHABLE | best probed damage 5.9057 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+` | UNREACHABLE | best probed damage 7.4894 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,ke+` | UNREACHABLE | best probed damage 3.2395 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `naive_E+` | UNREACHABLE | best probed damage 1.5360 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `naive_E+,ke+` | UNREACHABLE | best probed damage 4.6095 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |
| `naive_E-` | UNREACHABLE | best probed damage 3.1844 vs untreated 1.4575 — no improvement beyond the 14% noise floor at any potency |

## Within-dial order of the real arms, by the MRI channel (gap G4)

tau +0.80 over 10 within-dial pairs, permutation p = 0.024 (results/exam_v2.json S4); usable as a within-dial RANK input. Within-dial order of REAL arms by their trial's observed lesion ratio (backtest/potency.py). A rank input, not a magnitude: docs/RECOVERABILITY.md found the MRI-fitted potency biased 1.55x high. Never applied to survivors, which have no MRI trial.

| dial pattern | arms, lowest observed lesion ratio first |
|---|---|
| `ke` | ofatumumab (0.03, Gd-T1/scan), ocrelizumab (0.07, Gd-T1/scan) |
| `gamma_E` | natalizumab (0.17, new-T2), fingolimod (0.26, new-T2), cladribine (0.27, active-T2), ponesimod (0.44, CUAL/year) |
| `alpha_E` | dimethyl fumarate (0.15, new-T2), teriflunomide (0.20, Gd-T1/scan), IFN-beta (0.22, active-T2) |

Survivors are not in this table and cannot be: a screened candidate has no MRI trial. The order above is of drugs that already exist.

## Read this before quoting any line above

- **31 of 37 survivors claim a larger effect than the best drug ever tested in MS** (natalizumab, -68% relapse reduction in AFFIRM), and 22 of them claim better than -90%. Survivor effects here run -100% to -16%. **This is not credible and it is not meant to be read as a prediction.** The same model cannot reproduce a -30% effect for interferon beta on a real arm. The damage column measures the model, not the candidate.
- Nothing here is evidence about multiple sclerosis. validated=False throughout.
- Survivors are listed ALPHABETICALLY, never by predicted benefit. Ranking is gated on backtest/lomo.py beating its null; it does not (45.4pp vs 11.8pp).
- The model is blind to the depleting / sequestering / trafficking class — natalizumab, fingolimod, ponesimod, alemtuzumab and atacicept cannot come out beneficial in it at any potency (BUILD_PLAN §8.4). A candidate whose novelty is in trafficking is invisible to this screen.
- A surviving candidate has only avoided the failures this model can see.
