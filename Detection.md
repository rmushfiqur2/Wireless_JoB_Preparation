# Detection Theory: MAP/ML Detectors, Matched Filtering, Neyman-Pearson — Study Notes

Companion to the earlier notes — this section is the theoretical foundation *underneath* the BER analysis in the modulation/coding notes: how does a receiver actually decide what was sent, and in what sense is that decision "optimal"?

---

## 1. The Detection Problem, Formally

### 1.1 Setup

You observe `y`, generated under one of `M` hypotheses `H_0, H_1, ..., H_{M-1}` (e.g., which constellation symbol was sent). Each hypothesis specifies a conditional distribution `p(y | H_i)`. Detection theory asks: given `y`, which hypothesis should you declare, to optimize some criterion?

This is the exact theoretical machine underneath "pick the closest constellation point" from the modulation/coding notes — this section derives *why* that rule is optimal and under what assumptions.

### 1.2 Two different optimality criteria — the central distinction

| | **Maximum Likelihood (ML)** | **Maximum A Posteriori (MAP)** |
|---|---|---|
| Rule | `H_hat = argmax_i  p(y \| H_i)` | `H_hat = argmax_i  p(H_i \| y) = argmax_i  p(y \| H_i) * P(H_i)` |
| Uses prior `P(H_i)`? | No | Yes |
| Optimal when | Priors unknown, or equal priors | Priors known, minimizes **average** error probability |
| Relationship | MAP = ML when priors are equal | ML is a special case of MAP |

**The single most important fact to state cleanly:** MAP minimizes the *average* probability of error (over the prior distribution on hypotheses) — it is the **Bayes-optimal** detector under a 0/1 (uniform) cost/loss function. ML is what MAP reduces to when you either don't know the priors or the priors are equal — in that case the prior term `P(H_i)` is constant across hypotheses and drops out of the argmax, leaving pure likelihood comparison.

### 1.3 MAP derivation via Bayes' rule

```
p(H_i | y) = p(y | H_i) * P(H_i) / p(y)
```

Since `p(y)` doesn't depend on `i`, maximizing the posterior `p(H_i|y)` over `i` is equivalent to maximizing `p(y|H_i)*P(H_i)` over `i` — this is exactly the MAP rule. Be ready to write this out; it's a two-line derivation interviewers love to see done cleanly and confidently.

### 1.4 Binary case — the likelihood ratio test (LRT)

For `M=2` hypotheses, both ML and MAP reduce to a **likelihood ratio test**:

```
Lambda(y) = p(y|H_1) / p(y|H_0)   >< (threshold)
```

- **ML:** decide `H_1` if `Lambda(y) > 1`.
- **MAP:** decide `H_1` if `Lambda(y) > P(H_0)/P(H_1)` — the threshold shifts to account for the prior; if `H_0` is a priori more likely, you need *stronger* evidence (`y`) to switch your decision to `H_1`. This threshold-shifting intuition is a good one to be able to state in plain language, not just formula.

### 1.5 Worked example: binary antipodal signaling in AWGN (ties directly to the BPSK BER derivation from the modulation notes)

Hypotheses: `H_0: y = -s + n`, `H_1: y = +s + n`, with `n ~ N(0, sigma^2)`.

```
p(y|H_i) = (1/sqrt(2*pi*sigma^2)) * exp( -(y - mu_i)^2 / (2*sigma^2) )
```

Taking the log-likelihood ratio (monotonic transform, doesn't change the decision) and simplifying, the **ML/MAP-with-equal-priors decision rule collapses to:**

```
decide H_1 if y > 0,  decide H_0 if y < 0
```

i.e., **minimum-distance detection** — exactly the rule invoked (without derivation) in the modulation/coding notes. This is the clean, complete derivation that connects everything: the "pick the nearest constellation point" rule isn't an ad hoc heuristic, it *is* the ML/MAP detector for equally-likely, equal-energy signals in AWGN. Being able to state this chain of reasoning (Bayes rule -> LRT -> Gaussian likelihoods -> minimum distance) end-to-end is exactly the kind of connective, "sees the whole picture" answer a panel is listening for.

---

## 2. Matched Filtering

### 2.1 The problem matched filtering solves

Detection theory above assumed you already have a clean sample `y`. In practice, the receiver observes a continuous-time (or long discrete-time) signal corrupted by noise, and must first extract the best possible decision statistic from it. **Matched filtering** answers: *what linear filter, applied to the received waveform, maximizes SNR at the sampling instant, prior to making a detection decision?*

### 2.2 The matched filter result

Given a known transmitted pulse shape `s(t)` and AWGN, the filter that **maximizes output SNR** at the sampling instant is:

```
h(t) = s*(T - t)      (or in discrete time: h[n] = s*[N-1-n])
```

i.e., the filter impulse response is the **time-reversed, conjugated** copy of the transmitted pulse. Equivalently, matched filtering is **correlation** with the known signal shape:

```
y_out = integral y(t) * s*(t) dt      (correlator implementation, equivalent to the matched filter)
```

This is why matched filtering and **correlation receivers** are used interchangeably in textbooks — they're mathematically identical implementations of the same optimal operation.

### 2.3 Why it's optimal — the Cauchy-Schwarz argument (be ready to sketch this [Detailed Math](https://claude.ai/chat/572c4f50-0d35-48ea-8a41-ede3d1700fd6))

Output SNR for a linear filter `h(t)` applied to `s(t)+n(t)` (white noise, PSD `N0/2`) is maximized by Cauchy-Schwarz inequality applied to the correlation of `h` and `s` — equality (maximum) is achieved precisely when `h(t)` is proportional to `s*(T-t)`. The resulting **maximum achievable output SNR** is:

```
SNR_max = 2*E_s / N0
```

where `E_s = integral |s(t)|^2 dt` is the energy of the pulse. **Key insight:** the maximum output SNR depends *only* on the pulse's total energy, not its shape — this is a genuinely elegant and often-quoted result. It means matched filtering extracts *all* the available energy from the pulse regardless of its particular waveform shape, as long as you correlate against the exact known shape.

### 2.4 Matched filter and ML detection are the same operation

For AWGN channels, the matched filter output is a **sufficient statistic** for detection — running the matched filter and then doing minimum-distance/ML detection on the output is exactly equivalent to ML detection on the full continuous-time waveform. This is why real receivers are built as "matched filter (or equivalent correlator/RRC pulse-shaping filter) followed by a symbol-rate sampler and a minimum-distance detector" — that architecture *is* the ML receiver, not an approximation to it.

### 2.5 Practical instantiation: matched filtering in pulse-shaped systems

Real systems use **root-raised-cosine (RRC)** pulse shaping split between transmitter and receiver (each applies the square-root filter) so that transmit + receive filtering together form the full raised-cosine matched filter response, simultaneously controlling **inter-symbol interference (ISI)** — via the Nyquist zero-ISI criterion — and achieving the matched-filter SNR optimality. Worth mentioning if asked about *why* RRC filters specifically (rather than a simple rectangular matched filter) are used in practice: it's the matched-filter-plus-ISI-control combination.

---

## 3. Neyman-Pearson (NP) Detection

### 3.1 When MAP/ML aren't the right framework

MAP requires priors `P(H_i)`; ML implicitly assumes equal priors or doesn't care about them. But in some problems — **radar detection, spectrum sensing, anomaly/fault detection** — you genuinely don't have (or don't trust) a meaningful prior on `H_1` ("target present" / "signal present"), and more importantly, the *costs* of the two error types are wildly asymmetric (missing a real radar target vs. a false alarm have very different consequences). Neyman-Pearson handles exactly this case.

### 3.2 The NP criterion [Detailed Math](https://claude.ai/chat/93f20537-6620-41b4-aea8-825a0e3c6883)

Fix an acceptable **false alarm probability** `P_FA = Pr(decide H_1 | H_0 true) <= alpha`, and among all detectors satisfying that constraint, **maximize the probability of detection** `P_D = Pr(decide H_1 | H_1 true)`.

```
maximize P_D   subject to   P_FA <= alpha
```

### 3.3 The Neyman-Pearson Lemma (the key result [Detailed Math](https://claude.ai/chat/93f20537-6620-41b4-aea8-825a0e3c6883))

The optimal NP detector is **still a likelihood ratio test**:

```
decide H_1 if  Lambda(y) = p(y|H_1)/p(y|H_0)  >  eta
```

where the threshold `eta` is chosen specifically to satisfy `P_FA = alpha` exactly (not derived from priors/costs, but reverse-engineered from the false-alarm budget). This is a beautiful and very quotable result: **ML, MAP, and NP are all likelihood ratio tests** — they differ only in *how the threshold is set* (unity for ML, prior-ratio for MAP, false-alarm-constrained for NP). Being able to say this one sentence crisply is a strong signal of real understanding rather than memorized formulas.

### 3.4 ROC curves

The tradeoff between `P_D` and `P_FA` as the threshold `eta` varies is visualized as a **Receiver Operating Characteristic (ROC) curve** (`P_D` vs. `P_FA`). Key properties:
- Always passes through `(0,0)` and `(1,1)`.
- A better detector (or higher SNR) bows further toward the top-left corner (high `P_D` at low `P_FA`).
- The diagonal `P_D = P_FA` represents random guessing (no information) — any usable detector must lie above it.

This is directly relevant to **spectrum sensing** (a topic squarely in Qualcomm's wheelhouse, especially given Damnjanovic's spectrum-sharing/NR-U background) — deciding whether a channel is occupied is a textbook NP problem: false alarms waste spectrum opportunity, missed detections cause interference to the incumbent, and you typically have a hard regulatory constraint on `P_FA` (or on missed-detection probability) rather than a clean prior on channel occupancy.

---

## 4. Quick-reference formula sheet

```
ML rule:              H_hat = argmax_i p(y|H_i)
MAP rule:              H_hat = argmax_i p(y|H_i)*P(H_i)
Binary LRT (ML):        decide H_1 if p(y|H_1)/p(y|H_0) > 1
Binary LRT (MAP):        decide H_1 if p(y|H_1)/p(y|H_0) > P(H_0)/P(H_1)
Matched filter:          h(t) = s*(T - t)   [time-reversed, conjugated pulse]
Max output SNR:          SNR_max = 2*E_s / N0    (depends only on pulse energy)
NP rule:                decide H_1 if Lambda(y) > eta,  eta set so P_FA = alpha
```

---

## 5. Practice prompts

1. Derive the MAP decision rule from Bayes' rule and explain why it minimizes average error probability.
2. Show that for equal priors, MAP reduces exactly to ML.
3. Work through the binary antipodal-in-AWGN example and show the ML/MAP rule collapses to "decide H_1 if y > 0" — connect this explicitly to minimum-distance detection from the modulation notes.
4. State the matched filter impulse response and explain, at a high level, why time-reversal-plus-conjugation maximizes output SNR (Cauchy-Schwarz intuition is enough — full proof not required).
5. Explain why the maximum achievable output SNR from a matched filter depends only on pulse energy and not pulse shape — what's the intuition?
6. Explain the Neyman-Pearson criterion in one sentence, and explain why it's the right framework for spectrum sensing rather than MAP.
7. Sketch what a "good" vs. "bad" ROC curve looks like and explain what moving along a single ROC curve (fixed detector, varying threshold) physically means vs. moving between two different ROC curves (e.g., changing SNR).
8. Connect this to your own research: DOA estimation / array processing problems often involve a detection step (is there a source at angle theta, or just noise?) before an estimation step (what's the exact angle/amplitude?). How would you frame that detection sub-problem in ML/MAP/NP terms?
