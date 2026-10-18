# Propranolol enhances motivation for cognitive stability in tremor-dominant Parkinson’s disease

Analysis code for a double-blind, placebo-controlled, within-subject **crossover study** of the
β-adrenergic antagonist **propranolol** (40 mg) in **tremor-dominant Parkinson's disease**. The
study asks whether blocking β-adrenergic (noradrenergic) signalling changes working-memory
performance, the willingness to exert cognitive effort, and task-evoked pupil responses.

Patients completed, on placebo and on propranolol:

- a **colour-wheel working-memory task** crossing *effort level* (low vs high load) with
  *task-type* (resisting distraction vs flexibly updating memory), with simultaneous **pupillometry**;
- a **cognitive effort-discounting task** measuring willingness to repeat each kind of effort for reward.

## Key findings

- Propranolol **increased willingness to exert cognitive effort**, most clearly for the
  distractor-resistance (cognitive-stability) demand.
- It **modestly improved working-memory accuracy**, with a reaction-time slowing specific to the same condition.
- It **did not change task-evoked pupil size**, even though pupil tracked effort and task-type.
- Exploratorily, **higher working-memory pupil predicted lower subsequent effort choice**.

## What this repository shows

A complete mixed-model analysis pipeline in R:

- **(Generalized) linear mixed-effects models** (`lme4` / `lmerTest`) for accuracy, reaction time and choice;
- **Bayesian confirmation** of the key effects (`brms`), reported as posterior medians, credible
  intervals and probability of direction;
- **Estimated marginal means / simple-effect contrasts** (`emmeans`) with deviation coding so
  interaction models yield interpretable marginal effects;
- **Eye-tracking signal preprocessing**: trackloss screening, blink interpolation and smoothing
  of the raw pupil signal (`gazer` / `saccades`);
- **Pupillometry analysis**: trial-wise and binned time-course analyses of the task-evoked pupil response;
- **Publication tables and figures** (`sjPlot` / `kableExtra`, `ggplot2` / `ggrain`).

## Repository layout

```
code/
  drug_effects.Rmd            Main analysis: all tables and figures, top to bottom
  make_tab.R                  Renders a fitted model (or list of models) to a table image
  getmode.R                   Mode helper
  findBetas.R                 Wald z / p-values from a fitted mixed model
  preprocessing/
    pupil_preprocessing.Rmd   Eye-tracking: trackloss, blink interpolation, smoothing
    pupInterp.R               Interpolation / smoothing helper used by the above
```

**Pipeline order.** `preprocessing/pupil_preprocessing.Rmd` runs first and turns the raw
eye-tracking traces into the smoothed per-trial pupil data; `drug_effects.Rmd` then reads those
(plus the behavioural and clinical derivatives) and produces every table and figure.

The analysis notebook runs in sections that mirror the paper:

1. **Setup & data loading** — clinical/behavioural/pupil derivative files, outlier screening.
2. **Participant characteristics** (Table 1).
3. **Drug effects on performance** — accuracy and reaction time (Table 2).
4. **Drug effects on motivation** — effort-discounting choices (Table 3).
5. **Pupil drug effects** and **pupil tracking of effort/task-type** (Table 4; pupil ribbons).
6. **Pupil predicts subsequent effort choice** (Table 5).
7. **Bayesian confirmation** and **exploratory individual-difference** analyses (supplementary tables).

## Data availability

The behavioural, pupillometry and clinical (MDS-UPDRS, demographic) data are **sensitive patient
data and are not included** in this repository. They are available from the authors on reasonable
request, subject to the study's data-sharing agreement and ethics approval. To run the code, place
the derivative files in `Data/derivatives/` under the repository root (the notebook reads from there).

## Reproducing the analysis

Requires R (≥ 4.2). Install the packages the notebook lists on first run, plus `brms`, `lmerTest`
and `emmeans`:

```r
install.packages(c("tidyverse","lme4","lmerTest","emmeans","brms","sjPlot",
                   "kableExtra","ggrain","ggpubr","readxl","Rmisc","DescTools","webshot"))
# table images use webshot: webshot::install_phantomjs()
```

The preprocessing notebook additionally uses `data.table`, `zoo`, and the eye-tracking packages
`saccades` and `gazer` (the latter from GitHub: `remotes::install_github("dmirman/gazer")`).

Then knit `code/drug_effects.Rmd` from the repository root (or set the `PDPHARM_DIR`
environment variable to the repository root). Bayesian models are fit once and cached; set the
`run_brms` flag in the notebook to refit them.

## Reference

Papadopetraki et al. (in preparation). *Propranolol enhances motivation for cognitive
stability in tremor-dominant Parkinson’s disease.*

## Author

Danae Papadopetraki
