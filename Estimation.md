# Estimation Theory: MMSE, MLE, and the Cramer-Rao Bound — Study Notes

Companion to the detection theory notes — where detection asks "which discrete hypothesis," estimation asks "what continuous-valued parameter." This is the theoretical backbone for channel estimation, DOA estimation, and your own GS-SBL work, so lean into this section.

---

## 1. Framing: Detection vs. Estimation

| | Detection | Estimation |
|---|---|---|
| Unknown | Discrete hypothesis `H_i` (finite set) | Continuous parameter `theta` (real/complex-valued) |
| Goal | Pick the right hypothesis | Get `theta_hat` close to true `theta` |
| Optimality measure | Error probability | Mean-squared error (or other loss) |
| Analogous tool | MAP / ML / Neyman-Pearson | MMSE / MAP estimator / MLE |

The two are structurally parallel — MAP detection and MMSE/MAP estimation both come from Bayesian decision theory, just with different loss functions (0/1 loss for detection, squared-error loss for MMSE estimation).

---

## 2. Maximum Likelihood Estimation (MLE)

### 2.1 Definition

Given observed data `y` depending on an unknown deterministic (non-random) parameter `theta`, the MLE picks the value of `theta` that makes the observed data *most probable*:

```
theta_hat_ML = argmax_theta  p(y; theta)
```

Note the semicolon notation `p(y; theta)` (not `p(y|theta)`) — this signals `theta` is treated as a **fixed, unknown constant**, not a random variable with a prior. This is the key philosophical difference from MAP/Bayesian estimation (Section 3): **MLE is the "classical" (frequentist) approach**, no prior required or assumed.

### 2.2 Log-likelihood

In practice you maximize the **log-likelihood** `ln p(y;theta)` instead (monotonic transform, same maximizer, but sums instead of products for independent samples, and derivatives are cleaner):

```
theta_hat_ML = argmax_theta  ln p(y;theta)
```

Found by setting the derivative (the **score function**) to zero:

```
d/d(theta) [ ln p(y;theta) ] = 0
```

### 2.3 Worked example: MLE of a Gaussian mean (canonical, be ready to reproduce cold)

`y_1,...,y_N` i.i.d. `~ N(theta, sigma^2)`, `sigma^2` known.

```
ln p(y;theta) = -N/2 * ln(2*pi*sigma^2) - (1/(2*sigma^2)) * sum_i (y_i - theta)^2
```

Differentiate w.r.t. `theta`, set to zero:

```
(1/sigma^2) * sum_i (y_i - theta) = 0   =>   theta_hat_ML = (1/N) * sum_i y_i
```

**The MLE of the mean is just the sample average.** This is a good sanity-check example to keep in your back pocket — if a derivation you're doing under pressure doesn't reduce to something this intuitive in the Gaussian case, you likely made an error.

### 2.4 Key properties of MLE (know these, they get asked directly)

- **Asymptotically unbiased & efficient:** as `N -> infinity`, `theta_hat_ML` becomes unbiased and its variance approaches the Cramer-Rao bound (Section 4) — i.e., MLE is asymptotically optimal.
- **Consistency:** `theta_hat_ML -> theta_true` as `N -> infinity` (converges in probability).
- **Invariance property:** if `theta_hat` is the MLE of `theta`, then `g(theta_hat)` is the MLE of `g(theta)` for any function `g` — a genuinely useful and often-quoted shortcut (e.g., MLE of SNR follows directly from MLE of signal/noise power, no separate derivation needed).
- **Not necessarily unbiased for finite N** — a classic gotcha: the MLE of Gaussian *variance* (not mean) is biased for finite `N` (divides by `N` instead of `N-1`) — worth having as a concrete counterexample if someone probes "is MLE always unbiased?"

---

## 3. Bayesian Estimation: MMSE and MAP Estimators

### 3.1 The Bayesian shift

Unlike MLE, Bayesian estimation treats `theta` as a **random variable** with a known prior `p(theta)`, and uses the posterior `p(theta|y)` (via Bayes' rule, exactly as in the detection notes) to form an estimate.

### 3.2 MMSE estimator

**Definition:** the estimator minimizing mean-squared error `E[(theta - theta_hat)^2]`.

**Key result (very quotable, be ready to state without derivation):** the MMSE estimator is the **posterior mean** [Detailed Math](https://claude.ai/chat/1a74201d-06be-4c90-901a-03e68eccea59):

```
theta_hat_MMSE = E[theta | y] = integral theta * p(theta|y) d(theta)
```

**Why (brief sketch):** minimizing `E[(theta-theta_hat)^2 | y]` over `theta_hat` for fixed `y` is a basic calculus problem — differentiate w.r.t. `theta_hat`, set to zero, and you get `theta_hat = E[theta|y]` directly (the mean minimizes expected squared deviation — same fact as "the mean minimizes MSE" from basic statistics, just applied conditionally on `y`).

### 3.3 MAP estimator

**Definition:** picks the *mode* (peak) of the posterior rather than its mean [Detailed Math](https://claude.ai/chat/11736ec8-02b5-44a8-909f-8ee2d78d07ef):

```
theta_hat_MAP = argmax_theta  p(theta|y) = argmax_theta  p(y|theta)*p(theta)
```

Same Bayes'-rule structure as MAP *detection* from the previous notes — just now `theta` ranges over a continuum instead of a finite hypothesis set, and "argmax" is over a continuous variable (often solved via calculus, `d/d(theta)=0`, rather than by comparing a finite list).

### 3.4 MMSE vs. MAP — when do they differ, and why it matters

- **They coincide** when the posterior `p(theta|y)` is **symmetric and unimodal** (e.g., Gaussian posterior) — mean = mode in that case. This happens often in linear-Gaussian problems (Gaussian prior, Gaussian likelihood -> Gaussian posterior), which is why in many "nice" problems MMSE and MAP give the identical answer.
- **They differ** for skewed or multimodal posteriors — MAP can sit far from the posterior mean if the posterior has a long tail or multiple peaks. Be ready to say: **MMSE is generally the more principled criterion if you truly care about squared error**, but MAP is often cheaper to compute (pure optimization, no integration) and is what's implicitly used in many "regularized" estimation problems (see 3.5).

### 3.5 Bridge to your own research — MAP and sparse recovery

This is a direct, honest connection worth being ready to make explicitly: **sparse Bayesian learning (your GS-SBL work) is fundamentally a MAP/hierarchical-Bayesian estimation problem.** A Laplacian or sparsity-promoting prior on the parameter, combined with a Gaussian likelihood, yields a MAP estimator that is exactly the classic **L1-regularized least squares** (LASSO) problem — i.e., "sparse regularization" in compressed sensing is literally MAP estimation under a sparsity-inducing prior. If asked to connect estimation theory to your research, this is the clean, technically correct answer: your GS-SBL framework replaces a fixed sparsity prior with a hierarchical Gaussian prior (with hyperparameters learned from the data), which is a more flexible Bayesian estimation approach than plain MAP/LASSO, but still sits squarely in this Bayesian estimation framework.

---

## 4. The Cramer-Rao Bound (CRB)

### 4.1 What it is

The CRB is a **fundamental lower bound on the variance of any unbiased estimator**. It answers: no matter how clever your estimator design is, how precise can you *possibly* be, given the data model? This is the estimation-theory analog of the Shannon capacity limit — a hard theoretical ceiling on performance that concrete algorithms are measured against.

### 4.2 The formula (scalar case)

For an unbiased estimator `theta_hat` of a deterministic parameter `theta`, based on data with likelihood `p(y;theta)`:

```
Var(theta_hat) >= 1 / I(theta)
```

where `I(theta)` is the **Fisher information**:

```
I(theta) = E[ ( d/d(theta) ln p(y;theta) )^2 ]  =  -E[ d^2/d(theta)^2 ln p(y;theta) ]
```

The two forms (squared score vs. negative expected second derivative) are equal under standard regularity conditions — both are worth knowing since different derivations use different forms, and being able to switch between them signals real fluency rather than memorization.

### 4.3 Intuition for Fisher information

Fisher information measures **how sharply the log-likelihood is peaked** around the true parameter — a sharper peak (higher curvature, larger negative second derivative) means small changes in `theta` cause big changes in how probable the data looks, so the data is very informative about `theta`, so precise estimation is possible (low CRB, high `I(theta)`). A flat log-likelihood (low curvature) means the data barely distinguishes nearby values of `theta`, so no estimator can do well (high CRB, low `I(theta)`).

### 4.4 Worked example: CRB for the Gaussian mean (pairs directly with the MLE example above)

Using `ln p(y;theta) = -N/2*ln(2*pi*sigma^2) - (1/(2*sigma^2))*sum(y_i-theta)^2` from Section 2.3:

```
d/d(theta) ln p = (1/sigma^2) * sum(y_i - theta)
d^2/d(theta)^2 ln p = -N/sigma^2
```

```
I(theta) = N/sigma^2     =>     CRB = sigma^2 / N
```

**Check against the MLE:** `theta_hat_ML = sample mean`, and `Var(sample mean) = sigma^2/N` — **the sample mean's variance exactly equals the CRB.** This is the clean way to demonstrate "efficiency": the MLE achieves the Cramer-Rao bound exactly (not just asymptotically) in this particular Gaussian case. This worked example is extremely quotable in an interview — it shows the MLE, CRB, and the concept of an "efficient estimator" all in one clean, memorable calculation.

### 4.5 Efficient estimators

An estimator that **achieves** the CRB with equality is called an **efficient estimator**. Not every problem has one (the bound may not be achievable), but when one exists, it is automatically the MVUE (**Minimum Variance Unbiased Estimator**) — the best possible unbiased estimator in the mean-squared-error sense.

### 4.6 CRB for vector parameters (relevant to your array/DOA work)

For a vector parameter `theta` (e.g., multiple DOAs, or a channel impulse response with several taps), the CRB generalizes to a **matrix inequality**:

```
Cov(theta_hat) >= FIM^-1     (matrix, in the positive-semidefinite sense)
```

where `FIM` (Fisher Information Matrix) has entries `[FIM]_{ij} = -E[ d^2 ln p(y;theta) / (d theta_i * d theta_j) ]`. The diagonal entries of `FIM^-1` give the CRB on each individual parameter's variance. **This is precisely the tool used to derive fundamental limits on DOA estimation accuracy** (e.g., the classical CRB for DOA estimation with a uniform linear array, as a function of SNR, number of snapshots, and array geometry/aperture) — an extremely natural and honest bridge to your own research, more direct than almost any other topic in this whole prep set.

### 4.7 Why CRB matters practically

- Used to **benchmark real estimators** (e.g., is your DOA algorithm — MUSIC, ESPRIT, your GS-SBL approach — close to the CRB, or is there room for improvement?).
- Used in **system design** to predict best-case performance before building anything (e.g., "given this SNR and this many antennas/snapshots, what's the best possible DOA resolution we could ever achieve?").
- A very natural question to expect: *"How does the CRB for your problem scale with SNR / number of measurements / array size, and does your method approach it?"* — have a qualitative answer ready even if you don't have the exact scaling law memorized (typically CRB improves — decreases — with more snapshots/SNR/aperture, often close to linearly in SNR or inversely with N).

---

## 5. Quick-reference formula sheet

```
MLE:                       theta_hat_ML = argmax_theta  ln p(y;theta)
MMSE estimator:             theta_hat_MMSE = E[theta | y]      (posterior mean)
MAP estimator:              theta_hat_MAP = argmax_theta  p(y|theta)*p(theta)
Fisher information:         I(theta) = E[(d/d(theta) ln p)^2] = -E[d^2/d(theta)^2 ln p]
CRB (scalar):               Var(theta_hat) >= 1/I(theta),  for unbiased theta_hat
CRB (vector):                Cov(theta_hat) >= FIM^-1
Gaussian mean example:       theta_hat_ML = sample mean,  CRB = sigma^2/N  (MLE achieves CRB exactly)
```

---

## 6. Practice prompts

1. Derive the MLE of the mean of an i.i.d. Gaussian sample, and separately derive the Fisher information and CRB for the same problem — confirm the MLE achieves the bound.
2. State the difference between MLE and MAP estimation in one or two sentences — what does each assume about `theta`?
3. Explain why MMSE and MAP estimators coincide for a symmetric unimodal posterior but can diverge otherwise — sketch a simple example where they would differ.
4. Explain Fisher information in plain language: what does "high curvature of the log-likelihood" mean for how precisely you can estimate a parameter?
5. What does it mean for an estimator to be "efficient"? Is every parameter estimation problem guaranteed to have an efficient estimator?
6. Sketch how the CRB generalizes to a vector of parameters, and explain, at a high level, how this would apply to bounding DOA estimation accuracy for multiple sources.
7. Connect this to your own research: explain how your GS-SBL approach fits into the MAP/Bayesian estimation framework, and how the sparsity-promoting prior relates to L1-regularized (LASSO-style) MAP estimation.
