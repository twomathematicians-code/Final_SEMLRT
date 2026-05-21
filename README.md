# Small-Sample Corrections for Likelihood Ratio Tests in Structural Equation Modeling

A Monte Carlo simulation study comparing seven small-sample correction methods for the chi-square likelihood ratio test (LRT) in SEM, across three levels of model complexity and twelve sample sizes.

## Overview

In structural equation modeling (SEM), the standard chi-square test statistic ($$T_{ML}$$) often inflates Type I error rates at small sample sizes. This simulation evaluates how well seven correction methods recover the nominal alpha level (5%) under the null hypothesis (correctly specified models).

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
| 1 | **Uncorrected LRT** ($$T_{ML}$$) | Base test statistic: (n-1) x $$F_{ML}$$ | O(1/N) |
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
| $$F_{ML}$$ | Minimum of the ML discrepancy function |
| *S* | Sample covariance matrix (*p* x *p*) |
| $$\hat{\Sigma}$$| Model-implied covariance matrix (*p* x *p*) |
| $$\alpha$$ | Nominal significance level (0.05) |
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

- $$\hat{c}$$ = the Satorra-Bentler mean scaling factor (from above)
- $$\hat{a}$$ = an additional variance adjustment term derived from sixth-order moments

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

$$CI(\hat{\beta}) = \hat{\beta} \pm 1.96 \cdot SE(\hat{\beta})$$

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

$$\text{seed}(i) = 20240000 +\text{model}\_\text{offset} + i$$

where:

$$\text{model}\_\text{offset} = (\text{model}\_\text{rank} - 1) \times |\mathcal{N}| \times R$$

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

## Complete Architecture of the R Code

### High-Level Flow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                    SECTION 1: MODEL DEFINITIONS                      │
│                                                                     │
│  1. Define lavaan syntax for 3 models (Simple CFA, Complex CFA,     │
│     Full SEM)                                                      │
│  2. Build model registry (models list with name, syntax, desc)     │
│  3. Dummy-fit each model at N=2000 to extract df                   │
│                                                                     │
│  Result: model_info (name, syntax, df, n_indicators per model)     │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│              SECTION 2: SIMULATION ENGINE (7 chunks)                │
│                                                                     │
│  ┌─────────────────────────────────────────────────────────┐       │
│  │  CHUNK 2.1: BASE (Uncorrected ML)                        │       │
│  │  For each model × each N × each rep:                     │       │
│  │    1. set.seed(global_seed)                               │       │
│  │    2. simulateData(syntax, n)  →  n × p matrix           │       │
│  │    3. sem(syntax, data, "ML")  →  fit object              │       │
│  │    4. fitMeasures(fit, "chisq") → unc_chisq               │       │
│  │       fitMeasures(fit, "pvalue") → unc_p                 │       │
│  │    5. Store row: [Model, N, Rep, seed, chisq, p]          │       │
│  │  Save → thesis_sim_base_<model>.rds                       │       │
│  └─────────────────────────────────────────────────────────┘       │
│           │                                                         │
│           │  (reads base RDS to get exact N, Rep, seed)            │
│           ▼                                                         │
│  ┌─────────────────────────────────────────────────────────┐       │
│  │  CHUNK 2.2: SWAIN CORRECTION                             │       │
│  │  For each row in base:                                   │       │
│  │    1. set.seed(row$iter_seed)     ← SAME data as base    │       │
│  │    2. simulateData(syn, row$N)   ← identical data        │       │
│  │    3. sem(syn, data, "ML")       ← refit                 │       │
│  │    4. chisqSmallN(fit, "swain")  → T_Swain, P_Swain     │       │
│  │    5. Store row: [N, Rep, T_Swain, P_Swain]              │       │
│  │  Save → thesis_sim_swain_<model>.rds                     │       │
│  └─────────────────────────────────────────────────────────┘       │
│           │                                                         │
│           ▼  (same pattern for chunks 2.3 - 2.7)                   │
│  ┌─────────────────────────────────────────────────────────┐       │
│  │  CHUNK 2.3: BARTLETT  → chisqSmallN(fit, "bartlett")   │       │
│  │  CHUNK 2.4: YUAN 2015 → chisqSmallN(fit, "yuan.2015")  │       │
│  │  CHUNK 2.5: SATORRA-BENTLER → sem(..., estimator="MLM") │       │
│  │  CHUNK 2.6: YUAN-BENTLER    → sem(..., estimator="MLR") │       │
│  │  CHUNK 2.7: BOLLEN-STINE     → sem(..., test="bollen.stine",│     │
│  │                                   bootstrap=200)          │       │
│  │           (N <= 300 only — computationally expensive)     │       │
│  └─────────────────────────────────────────────────────────┘       │
│           │                                                         │
│           ▼                                                         │
│  ┌─────────────────────────────────────────────────────────┐       │
│  │  CHUNK 2.8: COMBINE ALL RESULTS                         │       │
│  │  Read 7 RDS files per model → LEFT JOIN on (N, Rep)     │       │
│  │  For each row, compute:                                  │       │
│  │    Reject_Unc   = (unc_p    < 0.05)  → TRUE/FALSE/NA    │       │
│  │    Reject_Swain = (P_Swain  < 0.05)                      │       │
│  │    Reject_Bart  = (P_Bart   < 0.05)                      │       │
│  │    Reject_Yuan  = (P_Yuan   < 0.05)                      │       │
│  │    Reject_SB    = (P_SB     < 0.05)                      │       │
│  │    Reject_YB    = (P_YB     < 0.05)                      │       │
│  │    Reject_BS    = (P_BS     < 0.05)                      │       │
│  │  Save → thesis_simulation_raw_results.rds                │       │
│  └─────────────────────────────────────────────────────────┘       │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│               SECTION 3: TYPE I ERROR COMPUTATION                   │
│  group_by(Model_Name, N) %>%                                       │
│    summarise(                                                        │
│      RR_Unc   = mean(Reject_Unc) * 100,    ← % rejected            │
│      RR_Swain = mean(Reject_Swain) * 100,  ← % rejected            │
│      ...                                                             │
│    )                                                                 │
│  Result: summary_table (36 rows = 3 models × 12 N values)          │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│               SECTIONS 4-6: VISUALIZATION & ANALYSIS                │
│  4.1 Line plots (log scale) with Bradley zone ribbon                │
│  4.2 Zoomed small-N plots                                            │
│  4.3 Bar chart comparison at key N                                   │
│  5.1-5.5 Convergence rate analysis (power-law fit on log-log)      │
│  6.1 Heatmap, 6.2 Bradley classification, 6.3 MAD ranking,         │
│  6.4 Interaction plot (df × method × N)                             │
└─────────────────────────────────────────────────────────────────────┘
```

---

### Detailed Per-Iteration Flow (Inside One Replication)

For one single replication (e.g., Simple CFA, N=100, Rep=42), the exact sequence is:

```
Step 1: SEED
    global_seed = 20240000 + (model_offset) + idx
    set.seed(global_seed)
    → This guarantees the SAME random data is generated every time
      this seed is used — critical for fair comparison across methods.

Step 2: DATA GENERATION
    sim_data <- simulateData(syntax, sample.nobs = n)
    |
    |  What this does internally:
    |  a) Parse the lavaan syntax string
    |  b) Construct implied covariance matrix Sigma_0 from the model
    |     (uses default starting values: loadings=1, residual variances=0.5, etc.)
    |  c) Draw n observations from N(0, Sigma_0) using Cholesky decomposition
    |  d) Return an n × p data.frame (e.g., 100 × 9 matrix for Simple CFA)
    |
    Output: sim_data = a 100×9 matrix of continuous normal data

Step 3: MODEL FITTING
    fit_ml <- sem(syntax, data = sim_data, estimator = "ML")
    |
    |  What this does internally:
    |  a) Compute sample covariance matrix S from sim_data (p×p)
    |  b) Optimize parameters (Lambda, Phi, Psi) to minimize:
    |        F_ML = log|Sigma_hat| + tr(S * Sigma_hat_inv) - log|S| - p
    |     where Sigma_hat = Lambda %*% Phi %*% t(Lambda) + Psi
    |  c) At convergence, compute test statistic:
    |        T_ML = (n - 1) * F_ML
    |  d) Compare T_ML to chi^2(df) → p-value
    |
    Output: fit_ml (lavaan object containing all estimates)

Step 4: EXTRACT STATISTICS
    unc_chisq <- fitMeasures(fit_ml, "chisq")   → e.g., 28.74
    unc_p     <- fitMeasures(fit_ml, "pvalue")  → e.g., 0.2321

Step 5: STORE
    results[[idx]] <- data.frame(
      Model = "simple_cfa",
      N = 100,
      Rep = 42,
      iter_seed = 20240042,
      unc_chisq = 28.74,
      unc_p = 0.2321
    )
```

Then for the Swain correction chunk, the **exact same seed is replayed**:

```
    set.seed(20240042)                       ← SAME seed
    sim_data <- simulateData(syn, 100)       ← IDENTICAL 100×9 data
    fit_ml <- sem(syn, sim_data, "ML")       ← IDENTICAL fit
    res <- chisqSmallN(fit_ml, "swain")      ← BUT different correction applied
    → T_Swain = T_ML / c_swain               ← divides by correction factor
    → P_Swain = pchisq(T_Swain, df)          ← new p-value
```

---

### The 7 Correction Methods — How Each One Works Internally

| # | Method | Code Mechanism | What It Does Mathematically |
|---|---|---|---|
| 1 | Uncorrected LRT | `fitMeasures(fit, "chisq")` | T_ML = (n-1) * F_ML. Compared directly to chi^2(df). Known to be inflated (too liberal) in small samples because the chi-square approximation is only asymptotic. |
| 2 | Swain | `chisqSmallN(fit, "swain")` | Applies Swain's correction factor c = 1 - (2t+2p+3)/(2n-2t-2p-1). Computes T_corrected = T_ML / c, then pchisq(T_corrected, df). c < 1, so T_corrected > T_ML, making it harder to reject (more conservative). |
| 3 | Bartlett (General) | `chisqSmallN(fit, "bartlett")` | Uses Bartlett's general formula derived from the trace of the information matrix. Similar to Swain but uses a different correction factor computation based on the full Edgeworth expansion. |
| 4 | Yuan (2015) | `chisqSmallN(fit, "yuan.2015")` | Yuan's third-order correction. Accounts for not just the first-order bias (like Bartlett) but also second-order terms from the Edgeworth expansion. Theoretically converges as O(1/n^3) vs O(1/n^2) for Bartlett. |
| 5 | Satorra-Bentler | `sem(..., estimator="MLM")` then `fitMeasures(fit, "chisq.scaled")` | Uses a mean-and-variance adjusted test statistic. Computes a scaling factor from fourth-order moments of the data. Designed for non-normality, but the scaling also helps in small samples. |
| 6 | Yuan-Bentler | `sem(..., estimator="MLR")` then `fitMeasures(fit, "chisq.scaled")` | Similar to Satorra-Bentler but uses a different (robust) scaling approach. MLR = Maximum Likelihood with Robust (Huber-White) standard errors. The scaled chi-square is T_ML / c where c accounts for both skewness and kurtosis. |
| 7 | Bollen-Stine | `sem(..., test="bollen.stine", bootstrap=200)` | Non-parametric bootstrap approach. For each of 200 bootstrap samples: transform data to fit the model, refit, collect T values. The p-value = proportion of bootstrap T values >= observed T. No distributional assumptions at all. Computationally expensive (~200x slower). |

---

### The Decision-Making Chain (Per Data Point)

For each of the 36,000 data points (3 models x 12 N x 1,000 reps), the code produces this decision chain:

```
                        ┌──────────────────────┐
                        │   p-value computed    │
                        │   for each method     │
                        └──────────┬───────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │    safe_reject(p, 0.05)    │
                    │                             │
                    │  if (length(p) != 1)        │──────► NA (convergence failure)
                    │    return NA                │
                    │  if (is.na(p))              │──────► NA (estimation failed)
                    │    return NA                │
                    │  if (p < 0.05)             │──────► TRUE  (REJECT H0)
                    │    return TRUE              │
                    │  else                       │──────► FALSE (FAIL TO REJECT H0)
                    │    return FALSE             │
                    └──────────────┬──────────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │   Aggregate over 1000 reps  │
                    │                             │
                    │  Rejection Rate =           │
                    │  mean(Reject == TRUE)       │
                    │  × 100  (in percent)        │
                    └──────────────┬──────────────┘
                                   │
                    ┌──────────────▼──────────────┐
                    │   CLASSIFY THE OUTCOME      │
                    │                             │
                    │  ┌─────────────────────┐    │
                    │  │ 4.5% ≤ RR ≤ 5.5%   │───► ROBUST
                    │  │ (Bradley's zone)    │    │
                    │  ├─────────────────────┤    │
                    │  │ RR < 4.5%           │───► CONSERVATIVE
                    │  │ (too few rejections) │    │
                    │  ├─────────────────────┤    │
                    │  │ RR > 5.5%           │───► LIBERAL
                    │  │ (too many rejections)│    │
                    │  └─────────────────────┘    │
                    │                             │
                    │  Ideal: RR = 5.0%           │
                    │  Tolerable: 4.5% - 5.5%     │
                    │  (Bradley, 1978)            │
                    └─────────────────────────────┘
```

---

### Complete Outcome Decision Matrix

| Outcome | RR Range | Meaning | Diagnosis | Action |
|---|---|---|---|---|
| **ROBUST** | 4.5% -- 5.5% | Test is well-calibrated | The correction successfully adjusts the chi-square approximation for this N and model complexity | Accept this method for this sample size |
| **SLIGHTLY LIBERAL** | 5.5% -- 7.0% | Rejects H0 slightly too often | Small residual bias; the correction helps but doesn't fully compensate | Method is usable but not ideal; consider larger N |
| **SEVERELY LIBERAL** | > 7.0% | Rejects H0 far too often | The correction is insufficient for this sample size / model complexity | Method should NOT be used at this N |
| **SLIGHTLY CONSERVATIVE** | 2.5% -- 4.5% | Rejects H0 slightly too rarely | Over-correction; the method is being too cautious | Acceptable for confirmatory research (low Type I error) |
| **SEVERELY CONSERVATIVE** | < 2.5% | Almost never rejects H0 | Gross over-correction; the correction is too aggressive | Method is problematic at this N |
| **NA (CONVERGENCE FAILURE)** | -- | Model didn't converge | N is too small relative to model complexity (n/p ratio too low) | Increase N or simplify model |

---

### Expected Pattern at Each N Level

```
N=50:    [████████████████████████████████] Uncorrected ~40% (SEVERELY LIBERAL)
         [███████████]                      Swain ~15%       (MODERATELY LIBERAL)
         [████████████████]                 Bartlett ~20%   (LIBERAL)
         [█████]                             Yuan ~8%        (SLIGHTLY LIBERAL)
         [███████████]                      SB ~12%          (MODERATELY LIBERAL)
         [██████]                            YB ~10%         (SLIGHTLY LIBERAL)
         [████]                              BS ~5%           (ROBUST) ← best
              │
N=200:   [█████████]                        Uncorrected ~10% (SLIGHTLY LIBERAL)
         [████]                              Swain ~5%       (ROBUST) ← good
         [█████]                             Bartlett ~6%    (ROBUST) ← good
         [████]                              Yuan ~5%        (ROBUST) ← good
              │
N=1000:  [█████]                             Uncorrected ~5.5% (ROBUST - barely)
         [█████]                             ALL methods ~5%  (ALL ROBUST)
              │
N=10000: [█████]                             ALL methods ~5.0% (ALL ROBUST)
```

---

### Data Flow Through RDS Files

The code uses **21 separate RDS files** as intermediate storage (7 methods x 3 models):

```
After Chunk 2.1 (Base):
  thesis_sim_base_simple_cfa.rds     → 12,000 rows × 6 cols  (12 N × 1000 reps)
  thesis_sim_base_complex_cfa.rds    → 12,000 rows × 6 cols
  thesis_sim_base_full_sem.rds       → 12,000 rows × 6 cols

After Chunk 2.2 (Swain):
  thesis_sim_swain_simple_cfa.rds    → 12,000 rows × 4 cols  (N, Rep, T_Swain, P_Swain)
  thesis_sim_swain_complex_cfa.rds   → 12,000 rows × 4 cols
  thesis_sim_swain_full_sem.rds      → 12,000 rows × 4 cols

  ... (same pattern for bart, yuan, sb, yb)

After Chunk 2.7 (Bollen-Stine, N<=300 only):
  thesis_sim_bs_simple_cfa.rds       →  4,000 rows × 3 cols  (4 N × 1000 reps)
  thesis_sim_bs_complex_cfa.rds      →  4,000 rows × 3 cols
  thesis_sim_bs_full_sem.rds         →  4,000 rows × 3 cols

After Chunk 2.8 (Combine):
  thesis_simulation_raw_results.rds → ~32,000 rows × 20 cols (3 models merged)
    Columns: Model, N, Rep, iter_seed, unc_chisq, unc_p,
             T_Swain, P_Swain, T_Bart, P_Bart, T_Yuan, P_Yuan,
             T_SB, P_SB, T_YB, P_YB, P_BS,
             Model_Name, Model_Desc, Model_DF,
             Reject_Unc, Reject_Swain, Reject_Bart, Reject_Yuan,
             Reject_SB, Reject_YB, Reject_BS
```

---

### Why the Seed Strategy Matters

The key design decision: every correction method chunk **re-generates the SAME data** by replaying the same seed stored in `base_df$iter_seed`. This means:

```
Base chunk:     seed=20240001 → data_A → ML fit → unc_p = 0.03
Swain chunk:    seed=20240001 → data_A → ML fit → swain_p = 0.07
Bartlett chunk: seed=20240001 → data_A → ML fit → bart_p = 0.08
Yuan chunk:     seed=20240001 → data_A → ML fit → yuan_p = 0.05
SB chunk:       seed=20240001 → data_A → MLM fit → sb_p = 0.06
YB chunk:       seed=20240001 → data_A → MLR fit → yb_p = 0.06
BS chunk:       seed=20240001 → data_A → BS fit → bs_p = 0.05

→ All 7 methods evaluated on IDENTICAL data → fair comparison
→ Paired comparison is possible (e.g., "Swain rejects but Bartlett doesn't")
```

Without this, each method would see different random data, adding noise and making paired comparisons invalid.

---

### Helper Functions

| Function | Purpose | Code |
|---|---|---|
| `safe_fit()` | Wraps `sem()` in `tryCatch`; returns `NULL` on convergence failure | `tryCatch(do.call(sem, args), error = function(e) NULL)` |
| `safe_fm()` | Extracts a fit measure; strips `lavaan.vector` S3 class; returns `NA` on failure | `as.numeric(unclass(fitMeasures(fit, measure)))` |
| `safe_smallN()` | Calls `chisqSmallN()`; returns `list(chisq, pvalue)` with NA fallbacks | `tryCatch(chisqSmallN(fit, method), error = function(e) NULL)` |
| `safe_reject()` | Binary rejection decision; handles NULL, length-0, NA gracefully | `if (is.na(p) \|\| length(p) != 1) return(NA); p < alpha` |

---

### Reproducibility (Section 7)

- `sessionInfo()` output for full environment capture
- **Global seed**: `set.seed(2024)`
- **Per-iteration seeding**: Deterministic seed formula ensures identical data generation across all correction method chunks
- **Cached chunks**: R Markdown `cache = TRUE` avoids redundant re-computation

## Reproducibility

See the detailed seed strategy explanation in the Architecture section above.

- **Global seed**: `set.seed(2024)`
- **Per-iteration seeding**: Deterministic seed formula ensures identical data generation across all correction method chunks (see "Why the Seed Strategy Matters")
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

### Complete Outcome Decision Matrix

See the detailed decision chain diagram in the Architecture section above.

| Outcome | RR Range | Meaning | Diagnosis | Action |
|---|---|---|---|---|
| **ROBUST** | 4.5% -- 5.5% | Test is well-calibrated | The correction successfully adjusts the chi-square approximation for this N and model complexity | Accept this method for this sample size |
| **SLIGHTLY LIBERAL** | 5.5% -- 7.0% | Rejects H0 slightly too often | Small residual bias; the correction helps but doesn't fully compensate | Method is usable but not ideal; consider larger N |
| **SEVERELY LIBERAL** | > 7.0% | Rejects H0 far too often | The correction is insufficient for this sample size / model complexity | Method should NOT be used at this N |
| **SLIGHTLY CONSERVATIVE** | 2.5% -- 4.5% | Rejects H0 slightly too rarely | Over-correction; the method is being too cautious | Acceptable for confirmatory research (low Type I error) |
| **SEVERELY CONSERVATIVE** | < 2.5% | Almost never rejects H0 | Gross over-correction; the correction is too aggressive | Method is problematic at this N |
| **NA (CONVERGENCE FAILURE)** | -- | Model didn't converge | N is too small relative to model complexity (n/p ratio too low) | Increase N or simplify model |

### Interpreting Rejection Rates (Bradley's Zones)

| Condition | Classification | Interpretation |
|---|---|---|
| RR in [2.5%, 7.5%] | Within Bradley's liberal zone | Acceptable for most purposes |
| RR in [4.5%, 5.5%] | Within Bradley's strict robust zone | Excellent Type I error control |
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

### Expected Pattern Across Sample Sizes

See the ASCII bar chart in the Architecture section above for a visual representation of how rejection rates are expected to change from N=50 to N=10,000 for each method.

## Computational Notes

- **Total iterations**: 3 models x 12 sample sizes x 1,000 reps = 36,000 per method
- **Bollen-Stine**: Restricted to N <= 300 (6 sample sizes x 1,000 reps x 3 models = 18,000 iterations with 200 bootstrap draws each)
- **Estimated runtime**: Several hours on a standard machine (depends on CPU; most expensive is Bollen-Stine)
- **Memory**: Per-model RDS files are ~1-2 MB each; combined results ~10-15 MB

## References
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. *Journal of Statistical Software*, 48(2), 1-36.
- Bradley, J. V. (1978). Robustness? *British Journal of Mathematical and Statistical Psychology*, 31(2), 144-152.
- Swain, A. J. (1975). A class of estimators for the mean of a multivariate normal distribution. *Unpublished doctoral dissertation*, University of Adelaide.
- Yuan, K.-H. (2015). Improved differential equation modeling of the mean and covariance structure. *Multivariate Behavioral Research*, 50(2), 164-184.
- Satorra, A., & Bentler, P. M. (1994). Corrections to test statistics and standard errors in covariance structure analysis. *Psychological Methods*, 1(4), 390-417.
- Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures in structural equation models. *Sociological Methodology*, 22, 111-135.

## License

This project is part of a master's thesis. Please cite appropriately if used in academic work.

