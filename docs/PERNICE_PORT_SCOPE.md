# Scoping the Pernice 2020 port — not started, specified

Written 2026-09-20. This file exists so the decision to port, or not to port, is
made against a specification rather than an impression. **Nothing here has been
implemented.**

    Pernice S, Follia L, Maglione A, Pennisi M, Pappalardo F, Novelli F,
    Clerico M, Beccuti M, Cordero F, Rolla S.
    "Computational modeling of the immune response in multiple sclerosis using
     epimod framework." BMC Bioinformatics 2020;21(Suppl 17):550.
    PMID 33308135, PMC7734848, doi:10.1186/s12859-020-03823-9.

## Licence — settled before anything else

| artifact | status | usable here? |
|---|---|---|
| the paper + Additional file 1 | **CC BY** (confirmed via Europe PMC `license: cc by`) | **yes** — transcribe, cite, tune nothing |
| `github.com/qBioTurin/Multiple-Sclerosis` (nets, R analysis) | **no licence file at all** → all rights reserved | **no** — do not vendor, copy or adapt |
| `github.com/qBioTurin/epimod` (framework) | GPL (>=2), per its DESCRIPTION | **no** — this repo is Apache-2.0 |

So the port is a transcription from the CC BY article, exactly as
`bricks/qsp_velez.py` was transcribed from Vélez de Mendizábal 2011. No code
from either repository enters this one.

## Why this model and not the other two

`bricks/qsp_velez.py` collapses depletion, sequestration and trafficking
blockade onto one dial, `gamma_E` — five of twelve quantified arms. Two other
candidates were assessed tonight and rejected: Martinez-Pasamar 2013 is the same
four ODEs, and Gazola 2025 has compartments but one blockade parameter on one
transmigration term, plus no Treg population at all.

Pernice has **two compartments coupled through an explicit BBB place**, keeps
**Teff/Treg cross-regulation**, and — the detail that decides it — carries
**separate BBB passage rates for the two cell types**:

    pPass_BBB_teff = 0.005        Teff_pass_BBB
    pPass_BBB_treg = 0.45         Treg_pass_BBB

A 90x difference, calibrated, with regulatory cells crossing far more readily
than effectors. A drug that blocks trafficking is then not the same object as a
drug that depletes, and a drug that blocks trafficking *asymmetrically* is
expressible for the first time in this repo.

## What a port would and would not separate

| currently collapsed | after the port | why |
|---|---|---|
| ocrelizumab, alemtuzumab (depletion) vs natalizumab (transit block) | **separated** | peripheral killing transitions vs `Teff_pass_BBB` |
| daclizumab's benefit vs its harm | **separated** | the paper gives `DACkillTeff` and `DACkillTreg` as two independent parameters — this repo's harm hypothesis as dials rather than a flag |
| fingolimod vs natalizumab | **still collapsed** | lymph node and blood vessel are ONE compartment, so there is no egress step for fingolimod to block |

## The shape of the work

26 places, 55 transitions (40 mass-action, 15 general), from which the authors
derive **a 26-ODE system with 20 calibrated parameters**. Two parameter sets,
healthy and MS, which **differ in only two values** — `pTeff_Activation`
(0.015 → 0.018) and `pTreg_Activation` (4e-04 → 7e-05).

Calibrated parameters (Additional file 1, Table S1), MS column where it differs:

    pTeff_Activation    0.015 / 0.018 (MS)   pTreg_Activation  4e-04 / 7e-05 (MS)
    pTreg_Dup           0.006                pTeff_Dup         0.04
    pTeff_KillsODC      6e-04                pTrkTe            0.02
    pTekA               6e-04                pPass_BBB_treg    0.45
    pPass_BBB_teff      0.005                pNKkillsTeff      0.01
    pNK_prod_IFNg       0.03                 pNK_prod_IL10     0.045
    pIL17_BBB           0.0115               pIL10_BBB         0.0765
    pRemyelinization    0.01                 pIL10Consuption   0.09
    pIL17Consuption     0.03                 pIFNgConsuption   0.05
    Cifn                20                   CIL10             10

Fixed parameters (Table S2): `FromTimoREG` 0.317, `FromTimoEFF` 0.296, `NKdup`
and `NKDegradation` and the four T-cell death/NLT transitions all 1/24 h^-1,
`Treg_prod_IL10` 0.05556, `Teff_prod_IL17` 0.00895, `Teff_prod_IFNg` 0.0466,
`DACDegradation` 0.001444057.

The general transitions are not mass action and must be transcribed as written;
the activation ones carry an IFN-gamma feedback coefficient in [0.5, 1.5]:

    f_TeffActivation_out = pTeff_Activation * x_RestingTeff_out * x_Antigen
                                            * (0.5 + exp(-x_IFNg_out / Cifn))
    f_TeffActivation_in  = pTeff_Activation * x_RestingTeff_in  * x_ODC_l
                                            * (0.5 + exp(-x_IFNg_in  / Cifn))
    f_MemActivation      = 0 before the 2nd antigen injection, then
                           2 * pTeff_Activation * x_Antigen * (0.5 + exp(...))
    f_NKentry            = 0 if x_NK_out > 30, else 0.267 * (30 - x_NK_out)

ODC carries five discrete myelination levels with an irreversible floor, so the
damage readout is a level count, not a continuous scalar — which is a **better**
fit for a lesion-count readout than the port currently in use.

## The bound that should be read before the expectations

Added 2026-09-20, after the parallel session measured a third failure mode and
`scripts/dial_ceiling.py` bounded the representation itself.

A drug in this repo is a set of multipliers on named dials, so **every arm
sharing a dial pattern gets the same prediction**. Real outcomes inside a dial
group are not the same:

    alpha_E   IFN-beta -30.0%, teriflunomide -31.5%, dimethyl fumarate -53.0%   spread 23.0pp
    gamma_E   natalizumab -68.0%, fingolimod -55.0%, alemtuzumab -55.0%,
              ponesimod -30.5%, cladribine -57.6%                               spread 37.5pp
    ke        ocrelizumab -46.0%, ofatumumab -50.0%                             spread  4.0pp

Predict each arm from its own dial group and see what error survives
(`scripts/dial_ceiling.py`, no simulation). **Which variant you use changes the
answer by a factor of five, so all three ship:**

| variant | MAE | null | headroom |
|---|---|---|---|
| in sample, all 12 arms | 6.6pp | 10.6pp | 4.0pp |
| singleton groups dropped (glatiramer and daclizumab are alone on their dials and fitted exactly for free) | 7.9pp | 10.5pp | 2.6pp |
| **out of sample, within group** — predict each arm from the OTHER arms in its group, which is what a model must do | **10.9pp** | **11.7pp** | **0.8pp** |

The last row is the one that counts, and it is scored the way every other
scorer in this repo is scored, null included. So:

- the representation is **not incapable**: a perfect dial-level model still
  beats the null. That is worth stating, because the obvious reading of the
  three failure modes is that it cannot, and that reading is wrong;
- but the prize is **0.8pp**, before any simulation or fitting error is charged
  against it. The measured LOMO is 45.9pp.

The in-sample 4.0pp was this repo's first estimate and it flattered the
representation twice over — by fitting each group's mean including the arm
being predicted, and by carrying two arms that are alone on their dials. The
parallel session caught both; it reproduced the model side exactly and the two
runs differed only on the null protocol, which is now stated in the script (an
out-of-sample model must be compared against an out-of-sample null, or it is
charged for information the null gets free).

**What that means for this port specifically.** A richer model does not need to
be better, it needs to land within about **one percentage point of perfect** to
clear a predict-the-mean null on this arm set. Pernice's extra compartments address
reachability — one of the three failures — and buy nothing against the other
two. Budget accordingly, and consider whether growing the arm set is the
cheaper lever, since the headroom is a property of the arms, not of the model.
Specifically **more arms per dial** — that is what the out-of-sample variant is
starved of, with ten arms across three multi-member groups.

## Honest expectations, stated before anyone starts

1. **A bigger model is not evidence it will score better.** The last time this
   repo replaced a toy with a grounded stack the direction gate went from 10/17
   to 5/16. Write the LOMO into the acceptance criteria *before* porting, not
   after.
2. **It does not touch the null problem.** Every out-of-sample scorer here loses
   to predict-the-mean, and a 26-state model fitted on 16 subjects' cytokine
   counts has more ways to lose, not fewer.
3. **It is stochastic by design** (SSA / tau-leaping). The 128-seed cohort and
   14% noise floor this repo already measured will apply, and the state count is
   6x larger — budget the compute before, not after.
4. **The calibration is not ours.** 8 MS patients and 8 healthy donors, cytokine
   counts at 18 hours. Inheriting it means inheriting its width.

## If it is done, the order

1. Transcribe states and mass-action transitions; pin Table S1 and S2 values in
   a module constant block with the table reference beside each one.
2. Transcribe the general transitions verbatim, including the IFN-gamma
   coefficient's [0.5, 1.5] range.
3. Reproduce a published figure, not a number — as `qsp_velez` reproduces
   Vélez's Figure 3. The HD-vs-MS contrast in their Figure 4 is the target.
4. Only then map arms onto dials, and only then re-run LOMO.
