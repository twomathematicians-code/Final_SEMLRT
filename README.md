# Small-Sample Corrections for Likelihood Ratio Tests in Structural Equation Modeling

A Monte Carlo simulation study comparing seven small-sample correction methods for the chi-square likelihood ratio test (LRT) in SEM, across three levels of model complexity and twelve sample sizes.

## Overview

In structural equation modeling (SEM), the standard chi-square test statistic (T_ML) often inflates Type I error rates at small sample sizes. This simulation evaluates how well seven correction methods recover the nominal alpha level (5%) under the null hypothesis (correctly specified models).

**Design:** 3 models x 7 correction methods x 12 sample sizes x 1,000 replications

## Research Questions

1. How does **model complexity** affect the performance of small-sample corrections?
2. Which correction method achieves Type I error rates closest to the nominal 5% across different sample sizes?
3. At what minimum sample size does each method fall within Bradley's (1978) robust zone?

## Models

| Model | Type | Latent Factors | Indicators (p) | Approx. df | Factor Correlations |
|---|---|---|---|---|---|
| Simple CFA | Confirmatory Factor Analysis | 3 | 9 | 24 | 1 |
| Complex CFA | Confirmatory Factor Analysis | 6 | 18 | ~120 | 8 |
| Full SEM | Full Structural Model | 3 | 9 | ~21 | 1 (+ structural paths) |

## Correction Methods

| # | Method | Implementation | Theoretical Error Order |
|---|---|---|---|
| 1 | **Uncorrected LRT** (T_ML) | Base test statistic: (n-1) x F_ML | O(1/N) |
| 2 | **Swain (Bartlett)** | `semTools::chisqSmallN(method = "swain")` | O(1/N^2) |
| 3 | **Bartlett (General)** | `semTools::chisqSmallN(method = "bartlett")` | O(1/N^2) |
| 4 | **Yuan (2015)** | `semTools::chisqSmallN(method = "yuan.2015")` | O(1/N^3) |
| 5 | **Satorra-Bentler** | `lavaan::sem(estimator = "MLM")` | O(1/N) |
| 6 | **Yuan-Bentler** | `lavaan::sem(estimator = "MLR")` | O(1/N) |
| 7 | **Bollen-Stine Bootstrap** | `lavaan::sem(test = "bollen.stine")` | O(1/B) |

## Mathematical Formulas

All mathematical formulas used in the simulation are listed below. Variables are defined as follows:

| Symbol | Definition |
|---|---|
| *n* | Sample size |
| *p* | Number of observed (manifest) variables |
| *t* | Number of latent variables |
| *df* | Model degrees of freedom |
| *F*_ML | Minimum of the ML discrepancy function |
| *S* | Sample covariance matrix (*p* x *p*) |
| Σ̂ | Model-implied covariance matrix (*p* x *p*) |
| α | Nominal significance level (0.05) |
| *B* | Number of bootstrap draws (B = 200 for Bollen-Stine) |

---

### 1. Base Test Statistic (Uncorrected LRT)

The standard ML chi-square test statistic is computed as:

$$T_{ML} = (n - 1) \; F_{ML}$$

where the ML discrepancy function is:

$$F_{ML} = \ln|\hat{\Sigma}| - \ln|S| + \text{tr}(S\hat{\Sigma}^{-1}) - p$$

Under the null hypothesis H₀, the asymptotic distribution is:

$$T_{ML} \xrightarrow{d} \chi^2(df) \quad \text{as } n \to \infty$$

The corresponding *p*-value is:

$$p_{unc} = P(\chi^2_{df} > T_{ML})$$

**Theoretical error order:** O(1/n)

---

### 2. Swain (1975) Correction

The Swain correction multiplies T_ML by a factor *c*_Swain derived from the model dimensionality. Two formulations exist:

**Swain c₁** (preferred when *n* is sufficiently large):

$$c_{Swain}^{(1)} = 1 - \frac{2t + 2p + 3}{2n - 2t - 2p - 1}$$

> **Condition:** Requires *n* > 2*t* + 2*p* + 2

**Swain c₂** (always valid for finite *n*):

$$c_{Swain}^{(2)} = \frac{n - \frac{2p + 5}{6} - \frac{2t + p}{3}}{n}$$

The corrected test statistic is:

$$T_{Swain} = c_{Swain} \times T_{ML}$$

The *p*-value is computed against the same chi-squared reference distribution:

$$p_{Swain} = P(\chi^2_{df} > T_{Swain})$$

**Implementation:** `semTools::chisqSmallN(fit, smallN.method = "swain")`

**Theoretical error order:** O(1/n²)

---

### 3. General Bartlett Correction

The Bartlett correction divides T_ML by a factor *c*_B that removes the first-order finite-sample bias in the chi-square approximation:

$$T_{Bart} = \frac{T_{ML}}{c_B}$$

The general Bartlett correction factor is derived from the Edgeworth expansion of T_ML:

$$c_B = 1 + \frac{1}{n} \left( \bar{a}_1 \right)$$

where *ā*₁ depends on model complexity and involves fourth-order moments of the observed variables. In `semTools`, *c*_B is computed internally based on the fitted model's information matrix and parameter count.

The expected value of the corrected statistic satisfies:

$$E(T_{Bart}) \approx df$$

The *p*-value:

$$p_{Bart} = P\left(\chi^2_{df} > \frac{T_{ML}}{c_B}\right)$$

**Implementation:** `semTools::chisqSmallN(fit, smallN.method = "bartlett")`

**Theoretical error order:** O(1/n²)

---

### 4. Yuan (2015) Correction

Yuan's (2015) correction removes the first **two** orders of finite-sample bias using a cubic correction based on the third-order Edgeworth expansion:

$$T_{Yuan} = \frac{T_{ML}}{c_{Yuan}}$$

where the correction factor is:

$$c_{Yuan} = 1 + \frac{\bar{a}_1}{n} + \frac{\bar{a}_2}{n^2}$$

The coefficients *ā*₁ and *ā*₂ are derived from cumulants of the sample covariance matrix and the model-implied parameter structure. The second-order term *ā*₂/n² is what distinguishes Yuan's correction from the standard Bartlett correction.

The *p*-value:

$$p_{Yuan} = P\left(\chi^2_{df} > \frac{T_{ML}}{c_{Yuan}}\right)$$

**Implementation:** `semTools::chisqSmallN(fit, smallN.method = "yuan.2015")`

**Theoretical error order:** O(1/n³)

---

### 5. Satorra-Bentler (1994) — MLM Estimator

The Satorra-Bentler correction uses a **mean-scaled** chi-square statistic that accounts for non-normality in the data. The scaling factor is computed from the weight matrix:

$$T_{SB} = \frac{T_{ML}}{\hat{c}}$$

where the scaling factor is:

$$\hat{c} = \frac{\text{tr}\left[(\mathbf{W}^{-1} - \mathbf{W}^{-1}\mathbf{P}\mathbf{W}^{-1})\hat{\mathbf{U}}\right]}{df}$$

with:

- **W** = the estimated weight matrix from the robust (ADF) estimator
- **P** = the product of the model Jacobian and its inverse
- **Û** = the estimated asymptotic covariance of the sufficient statistics

The *p*-value:

$$p_{SB} = P(\chi^2_{df} > T_{SB})$$

Under normality, ĉ → 1, so T_SB → T_ML (recovers the standard statistic).

**Implementation:** `lavaan::sem(syntax, data, estimator = "MLM")`, extracted via `fitMeasures(fit, c("chisq.scaled", "pvalue.scaled"))`

**Theoretical error order:** O(1/n) for mean-scaled version

---

### 6. Yuan-Bentler (1999) — MLR Estimator

The Yuan-Bentler correction extends Satorra-Bentler by additionally adjusting for the **variance** of the test statistic, producing a mean-and-variance adjusted statistic:

$$T_{YB} = \bar{T}_{SB} \quad \text{where} \quad \bar{T}_{SB} = \frac{T_{ML}}{\bar{d}}$$

The mean-and-variance adjustment factor is:

$$\bar{d} = \frac{\hat{c}}{1 + \frac{\hat{a} - 2\hat{c} + \hat{c}^2}{2 \cdot df \cdot \hat{c}^2}}$$

where:

- ĉ = the Satorra-Bentler mean scaling factor (from above)
- â = an additional variance adjustment term derived from sixth-order moments

The *p*-value uses a shifted chi-squared reference distribution:

$$p_{YB} = P(\chi^2_{\bar{df}} > \bar{d} \times T_{ML})$$

$$\bar{df} = \frac{df}{\bar{d}} \quad \text{(adjusted degrees of freedom)}$$

**Implementation:** `lavaan::sem(syntax, data, estimator = "MLR")`, extracted via `fitMeasures(fit, c("chisq.scaled", "pvalue.scaled"))`

**Theoretical error order:** O(1/n) for mean-and-variance adjusted version

---

### 7. Bollen-Stine (1992) Bootstrap

The Bollen-Stine bootstrap constructs an empirical reference distribution rather than relying on asymptotic theory:

**Step 1:** Transform the sample data so that it perfectly satisfies H₀:

$$\mathbf{S}^* = \hat{\Sigma} + \mathbf{S} - \hat{\Sigma} = \mathbf{S}$$

$$\mathbf{X}^{*}_{\text{pop}} \sim N(\mathbf{0}, \hat{\Sigma})$$

**Step 2:** For each bootstrap draw *b* = 1, ..., *B*:

$$\mathbf{X}^{(b)} \sim N(\mathbf{0}, \hat{\Sigma})$$

$$T^{(b)}_{ML} = \text{fit model to } \mathbf{X}^{(b)} \text{ and extract } T_{ML}^{(b)}$$

**Step 3:** Compute the bootstrap *p*-value:

$$p_{BS} = \frac{1}{B} \sum_{b=1}^{B} \mathbb{I}\left(T^{(b)}_{ML} \geq T_{ML}^{\text{obs}}\right)$$

where **I**(·) is the indicator function and *T*_ML^obs is the test statistic from the original sample.

**Implementation:** `lavaan::sem(syntax, data, test = "bollen.stine", bootstrap = 200)`, extracted via `fitMeasures(fit, "pvalue")`

**Theoretical error order:** O(1/B), where *B* = number of bootstrap draws

---

### 8. Type I Error Rate (Rejection Rate)

Since the data-generating model equals the fitted model (H₀ is true), the rejection rate estimates the Type I error rate. For each Model × *N* × Method combination:

$$\widehat{RR} = \frac{1}{R} \sum_{r=1}^{R} \mathbb{I}(p_r < \alpha)$$

where:

- *R* = number of replications (1,000)
- *p*_r = *p*-value from replication *r*
- α = 0.05 (nominal level)
- **I**(*p*_r < α) = 1 if *p*_r < α, 0 otherwise

The standard error of the estimated rejection rate:

$$SE(\widehat{RR}) = \sqrt{\frac{\widehat{RR}(1 - \widehat{RR})}{R}}$$

At *R* = 1,000 and α = 0.05: SE ≈ 0.69 percentage points.

---

### 9. Bradley's (1978) Robustness Zones

Bradley's criterion classifies a method's Type I error control quality:

| Zone | Bounds | Formula |
|---|---|---|
| **Liberal zone** | [α, 1.5α] | 0.05 ≤ RR ≤ 0.075 |
| **Robust zone** | [0.5α, 1.5α] | 0.025 ≤ RR ≤ 0.075 |
| **Strict robust zone** | [0.9α, 1.1α] | 0.045 ≤ RR ≤ 0.055 |

Classification rule used in the code:

$$\text{Class} = \begin{cases} \text{Conservative} & \text{if } RR < 0.025 \\ \text{Robust} & \text{if } 0.025 \leq RR \leq 0.075 \\ \text{Liberal} & \text{if } RR > 0.075 \end{cases}$$

> **Note:** The heatmap and interaction plots use the stricter [0.045, 0.055] band for visual emphasis.

---

### 10. Mean Absolute Deviation (MAD)

MAD measures the average absolute deviation of the rejection rate from the nominal 5% across all sample sizes:

$$\text{MAD}_{\text{all}} = \frac{1}{|\mathcal{N}|} \sum_{n \in \mathcal{N}} |RR(n) - 5|$$

where *N* is the set of all 12 sample sizes. Additional restricted MADs:

$$\text{MAD}_{N \leq 300} = \frac{1}{|\{n \in \mathcal{N} : n \leq 300\}|} \sum_{\substack{n \in \mathcal{N} \\ n \leq 300}} |RR(n) - 5|$$

$$\text{MAD}_{N=50} = |RR(50) - 5| \qquad \text{MAD}_{N=100} = |RR(100) - 5|$$

Methods are ranked by MAD (lower = better Type I error control).

---

### 11. Convergence Rate: Power-Law Regression

The empirical convergence rate estimates how fast each method's rejection rate approaches α as *n* → ∞. The power-law model:

$$\ln|RR(n) - \alpha \times 100| = \ln(a) - \beta \cdot \ln(n)$$

Equivalently:

$$|RR(n) - 5| = a \cdot n^{-\beta}$$

The estimated slope β̂ is obtained via OLS:

$$\hat{\beta} = -\frac{\sum_{i}(\ln n_i - \bar{n}_L)(\ln d_i - \bar{d}_L)}{\sum_{i}(\ln n_i - \bar{n}_L)^2}$$

where *d*_i = |RR(*n*_i) - 5| is the deviation at sample size *n*_i, and *n̄*_L and *d̄*_L are the log-space means.

Points with |RR - 5| ≤ 0.3% are excluded to avoid log-space instability near zero.

The 95% confidence interval:

$$CI_{95\%}(\hat{\beta}) = \hat{\beta} \pm 1.96 \cdot SE(\hat{\beta})$$

**Interpretation of β̂:**

| β̂ range | Convergence type | Corresponding order |
|---|---|---|
| β̂ < 1.5 | Linear | O(1/n) |
| 1.5 ≤ β̂ < 2.5 | Quadratic | O(1/n²) |
| β̂ ≥ 2.5 | Cubic or higher | O(1/n³) |

---

### 12. Tolerance Threshold: Minimum N\*

For each method, find the smallest sample size *n*\* where the rejection rate falls within tolerance ε of α:

$$N^* = \min\{n \in \mathcal{N} : |RR(n) - \alpha \times 100| \leq \varepsilon\}$$

Three tolerance levels are evaluated:

| ε | Interpretation |
|---|---|
| 0.25% | Very strict: RR must be in [4.75%, 5.25%] |
| 0.50% | Moderate: RR must be in [4.50%, 5.50%] |
| 1.00% | Lenient: RR must be in [4.00%, 6.00%] |

If no sample size achieves the tolerance, *N*\* is reported as NA.

---

### 13. Half-Decay Ratio

The half-decay ratio measures how quickly the deviation from 5% decreases between consecutive sample sizes:

$$\gamma_k = \frac{|RR(n_{k+1}) - 5|}{|RR(n_k) - 5|}$$

A value of γ_k < 0.5 means the deviation more than halved when moving from *n*_k to *n*_{k+1}. Smaller values indicate faster convergence.

---

### 14. Deterministic Seed Formula

To ensure identical data generation across all correction method chunks, each iteration uses a deterministic seed:

$$\text{seed}(i) = 20240000 + \text{model\_offset} + i$$

where:

$$\text{model\_offset} = (\text{model\_rank} - 1) \times |\mathcal{N}| \times R$$

- *i* = global iteration index (1 to |*N*| × *R* per model)
- model_rank = 1 (Simple CFA), 2 (Complex CFA), 3 (Full SEM)
- |*N*| = 12 sample sizes, *R* = 1,000 replications

This ensures that for iteration *i* of the Swain chunk, the same random data is generated as for iteration *i* of the Bartlett, Yuan, etc. chunks.

---

### 15. ML Discrepancy Function (Full Form)

For completeness, the full ML fitting function minimized by `lavaan::sem()` is:

$$F_{ML}(\theta) = \ln|\Sigma(\theta)| - \ln|S| + \text{tr}\left[S \cdot \Sigma(\theta)^{-1}\right] - p$$

where θ = vector of free model parameters. The parameter estimates are:

$$\hat{\theta} = \arg\min_{\theta} F_{ML}(\theta)$$

and the model-implied covariance matrix is:

$$\hat{\Sigma} = \Sigma(\hat{\theta})$$

---

### Formula Summary Table

| Formula | Section | Purpose |
|---|---|---|
| T_ML = (n−1) × F_ML | 1 | Base uncorrected test statistic |
| T_Swain = c_Swain × T_ML | 2 | Swain small-sample correction |
| T_Bart = T_ML / c_B | 3 | Bartlett correction |
| T_Yuan = T_ML / c_Yuan | 4 | Yuan (2015) higher-order correction |
| T_SB = T_ML / ĉ | 5 | Satorra-Bentler mean-scaled statistic |
| T_YB = T_ML / d̄ | 6 | Yuan-Bentler mean-and-variance adjusted |
| p_BS = (1/B) Σ I(T^(b) ≥ T_obs) | 7 | Bollen-Stine bootstrap p-value |
| RR = (1/R) Σ I(p_r < α) | 8 | Empirical Type I error rate |
| Bradley: 0.5α ≤ RR ≤ 1.5α | 9 | Robustness classification |
| MAD = mean(|RR − 5|) | 10 | Method ranking metric |
| ln|RR−5| = ln(a) − β ln(n) | 11 | Convergence rate estimation |
| N\* = min{n : |RR(n)−5| ≤ ε} | 12 | Minimum sample size threshold |
| γ_k = |dev(n_{k+1})| / |dev(n_k)| | 13 | Half-decay ratio |
| seed = 20240000 + offset + i | 14 | Reproducibility seed |

---

## Sample Sizes

```
50, 100, 150, 200, 250, 300, 500, 1000, 2000, 4000, 6000, 10000
```

- Bollen-Stine bootstrap is only run for N <= 300 (computational cost)
- Convergence validation is checked at N >= 1000

## Project Structure

```
.
├── thesis_simulation_lrt_corrections_4.1_final.Rmd   # Main analysis (knit to HTML)
├── README.md                                          # This file
└── (output RDS files from simulation chunks)           # Generated at runtime
```

### RDS Output Files (generated on knit)

| File | Description |
|---|---|
| `thesis_sim_base_*.rds` | Base uncorrected LRT results per model |
| `thesis_sim_swain_*.rds` | Swain correction results per model |
| `thesis_sim_bart_*.rds` | Bartlett correction results per model |
| `thesis_sim_yuan_*.rds` | Yuan (2015) correction results per model |
| `thesis_sim_sb_*.rds` | Satorra-Bentler (MLM) results per model |
| `thesis_sim_yb_*.rds` | Yuan-Bentler (MLR) results per model |
| `thesis_sim_bs_*.rds` | Bollen-Stine bootstrap results per model |
| `thesis_simulation_raw_results.rds` | Combined results across all methods |

## Code Structure & Execution Flow

### Phase 1: Setup

- Load packages: `lavaan`, `semTools`, `dplyr`, `tidyr`, `ggplot2`, `progress`
- Define global simulation parameters (sample sizes, replications, alpha)
- Define helper functions: `safe_fit()`, `safe_fm()`, `safe_smallN()`, `safe_reject()`

### Phase 2: Model Definitions (Section 1)

- Define lavaan syntax for three CFA/SEM models
- Build a model registry and compute degrees of freedom via dummy fits
- Store model metadata (name, description, syntax, df, number of indicators)

### Phase 3: Simulation Loop (Section 2)

Each correction method runs in a **separate cached R Markdown chunk** for fault isolation:

1. **Set seed** using a deterministic global index: `seed = 20240000 + model_offset + iter_index`
2. **Generate data** under H0 using `lavaan::simulateData()` with the correct model
3. **Fit model** using `lavaan::sem()` with the appropriate estimator/test
4. **Extract test statistic and p-value** via `fitMeasures()` or `chisqSmallN()`
5. **Store** results as a row in a data frame (N, Rep, seed, chi-sq, p-value)
6. **Save** per-model RDS files after each chunk completes

### Phase 4: Combine Results (Section 2.8)

- Merge all method-specific RDS files on `(N, Rep)`
- Compute binary rejection indicators for each method
- Save combined dataset as `thesis_simulation_raw_results.rds`

### Phase 5: Analysis (Sections 3-6)

1. **Rejection rates**: Proportion of p < 0.05 per Model x N x Method
2. **Convergence check**: Verify all methods approach 5% at N >= 1000
3. **Visualization**: Line plots, bar plots, focused small-N plots
4. **Convergence rate analysis**: Log-log power-law regression to estimate empirical convergence order
5. **Tolerance thresholds**: Minimum N* for accuracy at epsilon = 0.25%, 0.50%, 1.00%
6. **Bradley's robust zone classification**: Conservative / Robust / Liberal
7. **MAD ranking**: Mean Absolute Deviation from 5% for overall method comparison
8. **Interaction plots**: Method performance across model complexities

### Phase 6: Reproducibility (Section 7)

- `sessionInfo()` output for full environment capture

## Reproducibility

- **Global seed**: `set.seed(2024)`
- **Per-iteration seeding**: Deterministic seed formula ensures identical data generation across all correction method chunks
- **Cached chunks**: R Markdown `cache = TRUE` avoids redundant re-computation
- **Session info**: Full package versions recorded in Section 7

## How to Run

### Prerequisites

```r
install.packages(c("lavaan", "semTools", "dplyr", "tidyr",
                    "knitr", "progress", "ggplot2"))
```

### Run the Full Simulation

```r
# Knit the R Markdown document (produces HTML report)
rmarkdown::render("thesis_simulation_lrt_corrections_4.1_final.Rmd")
```

### Run Individual Sections

```r
# Or open in RStudio and run chunks interactively
# Each simulation chunk (2.1 - 2.7) is independent and cached
```

### Load Pre-computed Results

```r
results <- readRDS("thesis_simulation_raw_results.rds")
head(results)
```

## Key Output Tables & Plots

| Output | Section | Description |
|---|---|---|
| Rejection rate summary table | 3.1 | Full Model x N x Method rejection rates |
| Convergence check | 3.2 | Mean rejection rates at N >= 1000 |
| Main line plot | 4.1 | Rejection rates with Bradley's robust zone |
| Small-N focused plot | 4.2 | N <= 500 only |
| Bar comparison | 4.3 | Side-by-side bars at key sample sizes |
| Convergence rate table | 5.1 | Empirical beta estimates vs theoretical orders |
| Log-log convergence plot | 5.2 | Power-law fit visualization |
| Tolerance threshold table | 5.3 | Minimum N* for each method |
| Heatmap | 6.1 | Deviation from 5% (red = liberal, blue = conservative) |
| Bradley classification | 6.2 | Count of robust methods per condition |
| MAD ranking table | 6.3 | Method ranking by Mean Absolute Deviation |
| Interaction plot | 6.4 | Method x Model complexity interaction |

## Decision Framework

### Interpreting Rejection Rates

| Condition | Classification | Interpretation |
|---|---|---|
| RR in [2.5%, 7.5%] | Within Bradley's liberal zone | Acceptable for most purposes |
| RR in [4.5%, 5.5%] | Within Bradley's robust zone | Excellent Type I error control |
| RR < 2.5% | Severely conservative | Under-rejects; low statistical power |
| RR > 7.5% | Severely liberal | Over-rejects; inflated Type I error |

### Selecting a Correction Method

| Scenario | Recommended Approach |
|---|---|
| Very small N (< 100) | Bollen-Stine bootstrap or Yuan (2015) |
| Moderate N (100-300) | Swain or Bartlett correction |
| Large N (> 500) | Any method (all converge to 5%) |
| Non-normal data | Satorra-Bentler (MLM) or Yuan-Bentler (MLR) |
| High model complexity (df > 50) | Yuan (2015) for fastest convergence |
| Computational efficiency | Swain correction (closed-form, no iteration) |

## Computational Notes

- **Total iterations**: 3 models x 12 sample sizes x 1,000 reps = 36,000 per method
- **Bollen-Stine**: Restricted to N <= 300 (6 sample sizes x 1,000 reps x 3 models = 18,000 iterations with 200 bootstrap draws each)
- **Estimated runtime**: Several hours on a standard machine (depends on CPU; most expensive is Bollen-Stine)
- **Memory**: Per-model RDS files are ~1-2 MB each; combined results ~10-15 MB

## References

- Bradley, J. V. (1978). Robustness? *British Journal of Mathematical and Statistical Psychology*, 31(2), 144-152.
- Swain, A. J. (1975). A class of estimators for the mean of a multivariate normal distribution. *Unpublished doctoral dissertation*, University of Adelaide.
- Yuan, K.-H. (2015). Improved differential equation modeling of the mean and covariance structure. *Multivariate Behavioral Research*, 50(2), 164-184.
- Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. *Psychological Methods*, 1(4), 390-417.
- Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures in structural equation models. *Sociological Methodology*, 22, 111-135.
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. *Journal of Statistical Software*, 48(2), 1-36.

## License

This project is part of a master's thesis. Please cite appropriately if used in academic work.
