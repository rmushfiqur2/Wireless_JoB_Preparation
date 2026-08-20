# Common Distributions and Their Role in Wireless — Study Notes

Companion to the channel models and random processes notes — this section collects the specific distributions that keep reappearing throughout this whole prep series, with an explicit note on *where* each one already showed up and *why* it's the natural model for that situation.

---

## 1. Gaussian (Normal) Distribution

### 1.1 Definition

```
f(x) = 1/sqrt(2*pi*sigma^2) * exp( -(x-mu)^2 / (2*sigma^2) )
```

Parameterized by mean `mu`, variance `sigma^2`. Complex Gaussian: `Z = X+jY`, `X,Y` independent real Gaussians — this is the standard model for **thermal noise** (`n(t)` in every channel model in this series).

### 1.2 Why Gaussian shows up everywhere — the Central Limit Theorem (CLT)

The CLT is the single deepest reason Gaussian distributions dominate wireless modeling: **the sum of many independent, small random contributions tends toward Gaussian**, regardless of the individual contributions' own distributions. This single fact explains three separate distributions used across this prep series:

- **Thermal noise** is Gaussian because it results from the sum of enormous numbers of independent electron-thermal-motion contributions (physical CLT argument).
- **Shadowing** is log-normal (channel models notes, Section 3.1) because it's the *product* of many independent attenuation factors — a sum in the log/dB domain, hence Gaussian in dB, log-normal in linear scale.
- **Rayleigh fading** (Section 3 below) arises because the in-phase and quadrature components of a multipath sum are each (by the CLT) approximately Gaussian — Rayleigh is a *derived* distribution, not a primitive one, and its Gaussian origin is worth stating explicitly.

### 1.3 Where Gaussian already appeared

AWGN throughout the capacity/detection/estimation notes; MLE-of-the-mean worked example (estimation notes, Section 2.3); the "Gaussian input maximizes capacity" argument (capacity notes, Section 1.3); shadowing's log-normal derivation (channel models notes).

---

## 2. Exponential Distribution

### 2.1 Definition

```
f(x) = lambda * exp(-lambda*x),   x >= 0,   E[X] = 1/lambda
```

### 2.2 Role in wireless: instantaneous received power under Rayleigh fading

If the fading amplitude `|h|` is Rayleigh-distributed, the **instantaneous power/SNR** `|h|^2` is **exponentially distributed** — this is the single most important appearance of the exponential distribution in this whole prep series, since it's exactly what underlies the outage-capacity discussion (capacity notes, Section 2.3-2.4).

**Memoryless property:** the exponential is the *only* continuous distribution with the memoryless property, `Pr(X > s+t | X > s) = Pr(X > t)` — useful in queuing/traffic modeling (e.g., modeling call/session holding times, packet interarrival times in some traffic models) and worth knowing as a distinguishing mathematical fact if asked "why exponential and not some other distribution for this."

### 2.3 Closing the loop: "zero-outage capacity is zero"

Recall from the capacity notes: Rayleigh fading gives zero-outage capacity because the fading distribution has support down to zero with positive density. The exponential PDF `f(x) = lambda*exp(-lambda*x)` evaluated near `x=0` gives `f(0) = lambda > 0` — a nonzero density right at the origin — this is the precise, formal version of that earlier informal claim: **because the exponential density doesn't vanish at zero, there's always a nonzero probability the instantaneous SNR is arbitrarily close to zero, so no positive rate can be guaranteed with zero outage probability.**

---

## 3. Rayleigh and Rician Distributions

*(Full treatment in the channel models notes, Section 4.2-4.3 — summarized here for the "common distributions" list, with the derivation emphasis.)*

### 3.1 Rayleigh — the envelope of a zero-mean complex Gaussian

If `h = X+jY` with `X,Y` i.i.d. `N(0,sigma^2)` (no LOS component), the envelope `|h| = sqrt(X^2+Y^2)` is Rayleigh-distributed:

```
f(r) = (r/sigma^2)*exp(-r^2/(2*sigma^2)),   r >= 0
```

**Derivation intuition worth having ready:** this is a direct consequence of converting a 2D Gaussian (in Cartesian `X,Y` coordinates) into polar coordinates `(r,theta)` — the radius of a 2D isotropic Gaussian point is Rayleigh-distributed, and the phase `theta` comes out uniform on `[0,2*pi)`. This is a clean, self-contained derivation to be able to sketch if asked "where does Rayleigh actually come from."

### 3.2 Rician — envelope with a nonzero-mean (LOS) component

Same setup, but `X,Y` now have nonzero means (a deterministic LOS component plus Gaussian scatter) — the envelope follows a Rician distribution, governed by the K-factor as covered in the channel models notes.

---

## 4. Chi-Squared Distribution

### 4.1 Definition and connection to Gaussian

If `Z_1,...,Z_k` are i.i.d. standard Gaussian, `sum Z_i^2 ~ chi-squared with k degrees of freedom`. 

### 4.2 Role in wireless

- **`|h|^2` for Rayleigh fading is a chi-squared distribution with 2 degrees of freedom** (equivalently, an exponential — the two families overlap exactly at this special case: `chi-squared(2) = Exponential`). Worth knowing this equivalence explicitly, since it connects Sections 2 and 4 of this document into one fact rather than two.
- **Diversity combining:** with `L` independent diversity branches (e.g., `L`-branch spatial/antenna diversity or `L`-tap multipath combined via a RAKE receiver), the combined SNR after **maximum ratio combining (MRC)** follows a chi-squared distribution with `2L` degrees of freedom (sum of `L` independent exponential/chi-squared(2) variables) — this is the precise probabilistic statement behind "diversity improves the BER curve's slope" from the modulation/coding notes (Section 2.4): more diversity branches means the combined-SNR distribution concentrates more tightly around its mean (higher degrees of freedom -> less relative spread), directly producing the steeper BER-curve falloff associated with higher diversity order.
- **Detection theory connection:** chi-squared distributions appear naturally in energy-detector / spectrum-sensing test statistics (sum of squared Gaussian-noise samples under `H_0`) — directly relevant to the Neyman-Pearson spectrum-sensing discussion in the detection theory notes.

---

## 5. Nakagami-m Distribution (worth knowing as the "flexible" generalization)

### 5.1 Role

The **Nakagami-m distribution** is a flexible empirical fading model that generalizes both Rayleigh (`m=1` recovers Rayleigh exactly) and can approximate Rician fading for appropriate `m > 1`, while also being able to model **worse-than-Rayleigh fading** (`m < 1`) seen in some measured urban/indoor environments that Rayleigh doesn't fit well. It's popular in practice specifically because it has a convenient closed-form PDF that makes BER/outage analysis mathematically tractable across this whole range of fading severities, unlike Rician (whose PDF involves a Bessel function and is less analytically convenient). Worth mentioning as the practical/empirical counterpart to the "physically derived" Rayleigh/Rician pair — a good answer if asked "are there other fading models besides Rayleigh and Rician?"

---

## 6. Poisson Distribution and Poisson Process

### 6.1 Definition

```
P(X=k) = (lambda^k * exp(-lambda)) / k!
```

Discrete distribution for counting the number of events in a fixed interval, given events occur independently at a constant average rate `lambda`.

### 6.2 Role in wireless — two distinct uses worth distinguishing

- **Traffic modeling:** number of call/session/packet arrivals in a time interval — classical teletraffic theory (Erlang formulas for blocking probability in circuit-switched systems trace back to Poisson arrival assumptions).
- **Stochastic geometry / network topology modeling:** base stations or interferers/blockers scattered over an area are frequently modeled as a **Poisson Point Process (PPP)** — this is the modern tool used to derive tractable closed-form results for interference distribution, coverage probability, and (directly relevant to the channel models notes, Section 5.2) **mmWave blockage probability models**, where blockers (buildings, pedestrians) are modeled as a random spatial process and blockage probability is derived as a function of distance using PPP/Boolean-model machinery. Worth explicitly naming this as the mathematical tool underneath the "stochastic geometry approaches" mentioned informally in the channel models notes.

---

## 7. Uniform Distribution

### 7.1 Role in wireless

- **Carrier phase offset:** the phase of a received signal (relative to a reference) is often modeled as uniform on `[0,2*pi)` when there's no synchronization — this is exactly the phase term that comes out of the Rayleigh-derivation in Section 3.1 (polar decomposition of a 2D Gaussian yields uniform phase).
- **DOA priors in some Bayesian array-processing formulations:** without prior directional information, a uniform prior over the angular field of view is a natural non-informative choice — relevant if discussing Bayesian DOA methods (like GS-SBL) and what prior assumptions go into the model.

---

## 8. Summary Table — the "at a glance" version

| Distribution | Wireless role | Connects to |
|---|---|---|
| Gaussian | Thermal noise; CLT origin of others | AWGN capacity, MLE, matched filtering |
| Log-normal | Shadowing | Channel models notes Sec 3.1 |
| Rayleigh | Envelope of NLOS multipath fading | Channel models, outage capacity |
| Rician | Envelope of LOS + multipath fading | Channel models K-factor |
| Exponential | Instantaneous SNR/power under Rayleigh fading | Outage capacity, zero-outage-capacity result |
| Chi-squared | Diversity combining (MRC), energy detection | BER diversity slope, spectrum sensing |
| Nakagami-m | General empirical fading model | Generalizes Rayleigh/Rician |
| Poisson / PPP | Traffic arrivals; base station/blocker spatial modeling | mmWave blockage, interference modeling |
| Uniform | Carrier phase; non-informative DOA prior | Rayleigh derivation, Bayesian DOA |

---

## 9. Quick-reference formula sheet

```
Gaussian:            f(x) = 1/sqrt(2*pi*sigma^2) * exp(-(x-mu)^2/(2*sigma^2))
Exponential:          f(x) = lambda*exp(-lambda*x), x>=0,  E[X]=1/lambda
Rayleigh:            f(r) = (r/sigma^2)*exp(-r^2/(2*sigma^2))
Chi-squared(k):        sum of k squared i.i.d. standard Gaussians;  chi-sq(2) = Exponential
Poisson:             P(X=k) = lambda^k * exp(-lambda) / k!
Key equivalence:       Rayleigh envelope^2 = Exponential = Chi-squared(2 d.o.f.)
```

---

## 10. Practice prompts

1. Explain, via the Central Limit Theorem, why Gaussian, log-normal, and Rayleigh distributions are all "derived" consequences of the same underlying CLT argument, applied in slightly different ways.
2. Derive the fact that `|h|^2` is exponentially distributed when `h`'s envelope `|h|` is Rayleigh-distributed (simple change-of-variables exercise).
3. Explain how the memoryless property of the exponential distribution shows up (or doesn't) in fading-channel outage analysis.
4. Explain how diversity combining changes the distribution of combined SNR from chi-squared(2) to chi-squared(2L), and connect this to why higher diversity order produces a steeper BER-curve slope.
5. Explain why Nakagami-m is often preferred over Rician for tractable BER/outage analysis, despite Rician being more physically motivated for LOS scenarios.
6. Explain how a Poisson Point Process is used to model mmWave blockage probability, connecting back to the channel models notes.
7. Connect this to your own research: what distributional assumptions (e.g., on noise, on source amplitudes, on sparsity patterns) underlie your GS-SBL formulation, and how would you justify them if asked directly?
