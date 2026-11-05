# Alami et al. 2026 — Amazigh Nomads Show Strong Desire to Settle and Low Regret Despite Economic Costs

Replication code and manuscript source for the paper submitted to *Evolution and Human Sciences*. Analysis was run under R 4.6.1.

## Reproducing the paper

### Step 1 — Run the analysis

Open `Alami2026_ehs_analysis_code.R` in R and run it from top to bottom. This script:

- loads the data from `manuscript_data/`
- fits all statistical models (output saved to `model_output/`)
- generates all figures (saved to `figures/`)

The script must be run **with the repo root as the working directory** (i.e. `setwd()` to the folder containing this README, or open the repo as an RStudio project).

### Step 2 — Compile the manuscript

With [Quarto](https://quarto.org) installed, run from the repo root:

```bash
quarto render Alami2026_ehs.qmd
```

This produces `Alami2026_ehs.pdf`. Quarto requires a LaTeX distribution (e.g. [TinyTeX](https://yihui.org/tinytex/): `quarto install tinytex`).

The manuscript uses the `biblatex-chicago` LaTeX package (via the Cambridge journal extension). If the render fails with a missing-package error, install it via:

```bash
tlmgr install biblatex-chicago biber
```

Or from within R:

```r
tinytex::tlmgr_install(c("biblatex-chicago", "biber"))
```

---

## R package dependencies

Install all required packages before running the analysis script:

The Bayesian models (via `brms` and `cmdstanr`) require a C++ toolchain:

- **Windows**: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching your R version
- **Mac**: run `xcode-select --install` in Terminal
- **Linux**: install `r-base-dev` (included if R was installed via the CRAN apt repository)

```r
# cmdstanr is not on CRAN — install from Stan's r-universe
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
cmdstanr::install_cmdstan()

install.packages(c(
  "boot", "brms", "broom.helpers", "broom.mixed", "DiagrammeR", "dplyr", "ergm",
  "ggdist", "ggpattern", "ggplot2", "glmmTMB", "gt", "gtsummary",
  "janitor", "kableExtra", "kinship2", "lme4", "network", "patchwork",
  "scales", "sna", "statnet", "tibble", "tidyr",
  "bayestestR", "elevatr", "ggrepel", "ggridges",
  "ggspatial", "posterior", "rnaturalearth", "rnaturalearthdata",
  "sf", "terra", "tidyverse"
))
```

---

## Repository structure

```
.
├── Alami2026_ehs_analysis_code.R   # Main analysis script — run this first
├── Alami2026_ehs.qmd               # Manuscript source — compile this second
├── references.bib                  # Bibliography
├── variables.tex                   # LaTeX macros included by the QMD
├── code/                           # Helper scripts (see note below)
│   ├── dag_plot.R
│   ├── map_plot.R
│   ├── network_creation_function.R
│   ├── posterior_plots.R
│   └── Strand_models.R
├── manuscript_data/                # Input data (CSV)
├── figures/                        # Output figures (PNG) written by R script
├── model_output/                   # Model fit objects written by R script
└── _extensions/                    # Quarto Cambridge journal extension
```

The main R file sources only `dag_plot.R`, `map_plot.R` and `network_creation_function.R`. The manuscript's posterior and counterfactual figures are produced by the main R file itself, not by `posterior_plots.R`.

The income comparison figure (`figures/income_comparison.png`) uses a log(income + 1) axis labelled in US$, so that respondents with zero income stay in the plot.

---

## Data

The anonymised dataset is included in `manuscript_data/`. Raw interview records are not publicly available due to participant confidentiality.

## Funding

ENDOW (NSF awards 1743019, 2218860, and 2218861), Rep2SI (Leverhulme Award RL-2022-039), EIC (H2020 ERC grant #864519).
