# Screen results — candidates the model cannot rule out

> **RUN UNDER THE CAPACITY EXTENSION (K = 2000), NOT THE TRANSCRIPTION.**
> `bricks/qsp_velez.py`'s carrying capacity is an extension to Vélez de
> Mendizábal 2011, off by default. These results describe a DIFFERENT MODEL
> from `SCREEN_RESULTS.md`. The extension does not lift the ranking gate
> (`backtest/lomo_capacity.py`: 45.6pp against a 12.3pp null); it changes
> which candidates are screenable at all.

Generated 2026-09-21 by `PYTHONPATH=. python -m screen.report`. 98 candidates: every combination of up to 2 of the model's 7 named intervention points, each turned both ways at potency 0.5.

**This is a kill list with a remainder, not a ranking.** Survivors are listed ALPHABETICALLY, never by predicted benefit. Ranking is gated on backtest/lomo.py beating its null; it does not (45.4pp vs 11.8pp).

## Counts

| verdict | n |
|---|---|
| DEGENERATE | 8 |
| SURVIVED | 30 |
| UNREACHABLE | 60 |

## Survivors (30)

Alphabetical. The damage column is a simulated median over 8 infection histories, shown so a reader can challenge the verdict — not as a score.

| candidate | best simulated damage | vs untreated | at potency |
|---|---|---|---|
| `alpha_E+,gamma_R-` | 0.0018 | -97% | 0.8 |
| `alpha_E+,ke-` | 0.0088 | -84% | 0.8 |
| `alpha_E-,delta-` | 0.0390 | -28% | 0.8 |
| `alpha_E-,gamma_E+` | 0.0391 | -28% | 0.8 |
| `alpha_E-,gamma_R-` | 0.0014 | -97% | 0.8 |
| `alpha_E-,ke-` | 0.0052 | -90% | 0.8 |
| `alpha_E-,naive_E-` | 0.0425 | -21% | 0.8 |
| `alpha_R+,gamma_R-` | 0.0015 | -97% | 0.8 |
| `alpha_R+,ke-` | 0.0022 | -96% | 0.8 |
| `alpha_R+,naive_E+` | 0.0325 | -40% | 0.8 |
| `alpha_R+,naive_E-` | 0.0324 | -40% | 0.8 |
| `alpha_R-,gamma_R-` | 0.0025 | -95% | 0.8 |
| `delta+,ke-` | 0.0081 | -85% | 0.8 |
| `delta-,ke-` | 0.0033 | -94% | 0.8 |
| `gamma_E+,alpha_R+` | 0.0312 | -42% | 0.8 |
| `gamma_E+,delta-` | 0.0422 | -22% | 0.8 |
| `gamma_E+,gamma_R-` | 0.0007 | -99% | 0.8 |
| `gamma_E+,ke-` | 0.0055 | -90% | 0.8 |
| `gamma_E-,alpha_R+` | 0.0393 | -27% | as enumerated |
| `gamma_E-,gamma_R-` | 0.0225 | -58% | 0.8 |
| `gamma_E-,ke-` | 0.0237 | -56% | 0.8 |
| `gamma_R-` | 0.0016 | -97% | 0.8 |
| `gamma_R-,delta+` | 0.0016 | -97% | 0.8 |
| `gamma_R-,delta-` | 0.0015 | -97% | 0.8 |
| `gamma_R-,ke+` | 0.0023 | -96% | 0.8 |
| `gamma_R-,ke-` | 0.0013 | -98% | 0.8 |
| `gamma_R-,naive_E+` | 0.0040 | -93% | 0.8 |
| `gamma_R-,naive_E-` | 0.0004 | -99% | 0.8 |
| `naive_E+,ke-` | 0.0066 | -88% | 0.8 |
| `naive_E-,ke-` | 0.0061 | -89% | 0.8 |

## Killed

| candidate | verdict | why |
|---|---|---|
| `alpha_E+,alpha_R+` | DEGENERATE | same intervention points as lenercept |
| `alpha_E-,alpha_R+` | DEGENERATE | same intervention points as lenercept |
| `alpha_R+` | DEGENERATE | same intervention points as daclizumab |
| `alpha_R+,delta+` | DEGENERATE | same intervention points as glatiramer acetate |
| `alpha_R+,delta-` | DEGENERATE | same intervention points as glatiramer acetate |
| `delta-` | DEGENERATE | same intervention points as abatacept |
| `delta-,naive_E-` | DEGENERATE | same intervention points as IFN-gamma |
| `ke-` | DEGENERATE | same intervention points as rituximab |
| `alpha_E+` | UNREACHABLE | best probed damage 0.0552 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,alpha_R-` | UNREACHABLE | best probed damage 0.0697 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,delta+` | UNREACHABLE | best probed damage 0.0556 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,delta-` | UNREACHABLE | best probed damage 0.0485 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,gamma_E+` | UNREACHABLE | best probed damage 0.0549 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,gamma_E-` | UNREACHABLE | best probed damage 0.0548 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,gamma_R+` | UNREACHABLE | best probed damage 0.0744 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,ke+` | UNREACHABLE | best probed damage 0.0723 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,naive_E+` | UNREACHABLE | best probed damage 0.0550 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E+,naive_E-` | UNREACHABLE | best probed damage 0.0553 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-` | UNREACHABLE | best probed damage 0.0471 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,alpha_R-` | UNREACHABLE | best probed damage 0.0677 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,delta+` | UNREACHABLE | best probed damage 0.0476 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,gamma_E-` | UNREACHABLE | best probed damage 0.0530 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,gamma_R+` | UNREACHABLE | best probed damage 0.0720 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,ke+` | UNREACHABLE | best probed damage 0.0693 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_E-,naive_E+` | UNREACHABLE | best probed damage 0.0499 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R+,gamma_R+` | UNREACHABLE | best probed damage 0.0589 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R+,ke+` | UNREACHABLE | best probed damage 0.0613 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-` | UNREACHABLE | best probed damage 0.0687 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,delta+` | UNREACHABLE | best probed damage 0.0691 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,delta-` | UNREACHABLE | best probed damage 0.0681 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,gamma_R+` | UNREACHABLE | best probed damage 0.0964 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,ke+` | UNREACHABLE | best probed damage 0.0872 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,ke-` | UNREACHABLE | best probed damage 0.0527 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,naive_E+` | UNREACHABLE | best probed damage 0.0691 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `alpha_R-,naive_E-` | UNREACHABLE | best probed damage 0.0685 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta+` | UNREACHABLE | best probed damage 0.0545 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta+,ke+` | UNREACHABLE | best probed damage 0.0715 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta+,naive_E+` | UNREACHABLE | best probed damage 0.0544 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta+,naive_E-` | UNREACHABLE | best probed damage 0.0547 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta-,ke+` | UNREACHABLE | best probed damage 0.0702 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `delta-,naive_E+` | UNREACHABLE | best probed damage 0.0473 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+` | UNREACHABLE | best probed damage 0.0522 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,alpha_R-` | UNREACHABLE | best probed damage 0.0677 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,delta+` | UNREACHABLE | best probed damage 0.0537 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,gamma_R+` | UNREACHABLE | best probed damage 0.0728 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,ke+` | UNREACHABLE | best probed damage 0.0698 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,naive_E+` | UNREACHABLE | best probed damage 0.0528 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E+,naive_E-` | UNREACHABLE | best probed damage 0.0515 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-` | UNREACHABLE | best probed damage 0.0539 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,alpha_R-` | UNREACHABLE | best probed damage 0.0709 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,delta+` | UNREACHABLE | best probed damage 0.0543 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,delta-` | UNREACHABLE | best probed damage 0.0533 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,gamma_R+` | UNREACHABLE | best probed damage 0.0734 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,ke+` | UNREACHABLE | best probed damage 0.0718 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,naive_E+` | UNREACHABLE | best probed damage 0.0542 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_E-,naive_E-` | UNREACHABLE | best probed damage 0.0538 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+` | UNREACHABLE | best probed damage 0.0734 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,delta+` | UNREACHABLE | best probed damage 0.0740 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,delta-` | UNREACHABLE | best probed damage 0.0724 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,ke+` | UNREACHABLE | best probed damage 0.0923 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,ke-` | UNREACHABLE | best probed damage 0.0553 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,naive_E+` | UNREACHABLE | best probed damage 0.0736 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `gamma_R+,naive_E-` | UNREACHABLE | best probed damage 0.0732 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `ke+` | UNREACHABLE | best probed damage 0.0710 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `naive_E+` | UNREACHABLE | best probed damage 0.0540 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `naive_E+,ke+` | UNREACHABLE | best probed damage 0.0712 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `naive_E-` | UNREACHABLE | best probed damage 0.0541 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |
| `naive_E-,ke+` | UNREACHABLE | best probed damage 0.0707 vs untreated 0.0541 — no improvement beyond the 14% noise floor at any potency |

## Read this before quoting any line above

- **20 of 30 survivors claim a larger effect than the best drug ever tested in MS** (natalizumab, -68% relapse reduction in AFFIRM), and 15 of them claim better than -90%. Survivor effects here run -99% to -21%. **This is not credible and it is not meant to be read as a prediction.** The same model cannot reproduce a -30% effect for interferon beta on a real arm. The damage column measures the model, not the candidate.
- Nothing here is evidence about multiple sclerosis. validated=False throughout.
- Survivors are listed ALPHABETICALLY, never by predicted benefit. Ranking is gated on backtest/lomo.py beating its null; it does not (45.4pp vs 11.8pp).
- The model is blind to the depleting / sequestering / trafficking class — natalizumab, fingolimod, ponesimod, alemtuzumab and atacicept cannot come out beneficial in it at any potency (BUILD_PLAN §8.4). A candidate whose novelty is in trafficking is invisible to this screen.
- A surviving candidate has only avoided the failures this model can see.
