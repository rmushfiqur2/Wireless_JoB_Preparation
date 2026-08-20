# Random Processes: Stationarity, Autocorrelation, PSD — Study Notes

Companion to the channel models and detection theory notes — this is the mathematical machinery that made statements like "white noise," "coherence time," and "coherence bandwidth" precise in earlier sections. Framed here around where each concept already showed up, since this material is rarely tested in isolation.

---

## 1. What a Random Process Is

### 1.1 Definition

A **random process** `X(t)` (continuous-time) or `X[n]` (discrete-time) is a collection of random variables indexed by time — for each fixed time `t`, `X(t)` is a random variable; for each fixed outcome (sample path/realization), `X(t)` is a deterministic function of time. Noise `n(t)` and fading channel gain `h(t)` throughout the earlier notes are both random processes in exactly this sense.

### 1.2 Statistical description

Fully describing a random process in general requires the **joint distribution** of `X(t_1), X(t_2), ..., X(t_k)` for every possible set of times and every `k` — intractable in general. In practice, essentially everything in communications theory works with a much smaller set of summary statistics:

```
mean:            mu_X(t) = E[X(t)]
autocorrelation: R_X(t1,t2) = E[X(t1)*X*(t2)]
autocovariance:  C_X(t1,t2) = E[(X(t1)-mu_X(t1))*(X(t2)-mu_X(t2))*]
```

---

## 2. Stationarity

### 2.1 Strict-sense stationarity (SSS)

A process is **strictly stationary** if its complete joint statistics are invariant to a time shift — the joint distribution of `X(t_1),...,X(t_k)` equals that of `X(t_1+tau),...,X(t_k+tau)` for any shift `tau`. This is a strong condition, rarely verified directly in practice.

### 2.2 Wide-sense stationarity (WSS) — the one that actually gets used everywhere

A process is **wide-sense stationary** if only its first two moments are shift-invariant:

```
1. mu_X(t) = mu_X   (constant mean, doesn't depend on t)
2. R_X(t1,t2) = R_X(t1-t2) = R_X(tau)    (autocorrelation depends only on the time DIFFERENCE tau = t1-t2, not absolute time)
```

**Why WSS is the workhorse assumption:** almost every result in this whole prep series — matched filtering, capacity formulas, PSD-based bandwidth arguments — implicitly assumes WSS. It's a much weaker (easier to satisfy/verify) condition than SSS, but sufficient for everything built on autocorrelation/PSD, which is why you'll see "WSS" as a standing assumption rather than SSS almost everywhere in the comms literature.

**Key relationship:** SSS implies WSS (a stronger condition implies the weaker one), but not vice versa in general — the one common exception worth knowing: for a **Gaussian process**, WSS *does* imply SSS, because a Gaussian process's complete statistics are fully determined by its mean and covariance alone — no higher-order moments carry extra information. This is a genuinely elegant fact and a good one to have ready if asked "when does WSS give you the full picture?"

### 2.3 Where WSS already appeared (implicitly) in this prep series

- **AWGN in the capacity/detection notes:** "white" noise with PSD `N0/2` is a WSS assumption by construction — white noise's autocorrelation is a delta function (Section 3.2), which only makes sense to write down under a shift-invariance assumption.
- **Coherence time/bandwidth (channel models notes, Section 4.4-4.5):** these concepts are defined via the autocorrelation of the *fading process* `h(t)` — coherence time is literally "the `tau` at which the channel's autocorrelation function drops significantly," which presumes WSS (autocorrelation as a function of lag `tau` alone, not of absolute time).

---

## 3. Autocorrelation Function

### 3.1 Definition (WSS case)

```
R_X(tau) = E[ X(t+tau) * X*(t) ]
```

Measures how similar the process is to a time-shifted version of itself. Properties worth knowing cold:

- `R_X(0) = E[|X(t)|^2]` = **average power** of the process — a very commonly used fact (autocorrelation at zero lag *is* signal power).
- `R_X(tau)` is Hermitian symmetric: `R_X(-tau) = R_X*(tau)` (real and even, `R_X(-tau)=R_X(tau)`, for real-valued processes).
- `|R_X(tau)| <= R_X(0)` — maximum always at zero lag (the process is always at least as correlated with itself at zero shift as at any other shift).
- **Decorrelation time:** the lag `tau` beyond which `R_X(tau)` becomes negligible — this is exactly the formal definition underlying "coherence time" from the channel models notes: coherence time is the decorrelation time of the fading process.

### 3.2 White noise — the canonical extreme case

"White" noise is defined by an autocorrelation that is a **delta function**:

```
R_N(tau) = (N0/2) * delta(tau)
```

**Physical meaning:** the noise is *completely* uncorrelated with itself at any nonzero time lag, no matter how small — every sample is statistically independent of every other sample (for Gaussian noise, uncorrelated implies independent). This is the precise mathematical statement behind "AWGN" used throughout the capacity, detection, and estimation notes — worth being able to state that "white" specifically refers to this delta-function autocorrelation property, not just "Gaussian" (Gaussian describes the *marginal distribution*; white describes the *correlation structure* — they're independent properties, and AWGN combines both).

---

## 4. Power Spectral Density (PSD)

### 4.1 The Wiener-Khinchin theorem

The PSD is the Fourier transform of the autocorrelation function:

```
S_X(f) = integral R_X(tau) * exp(-j*2*pi*f*tau) d(tau)
```

**This single relationship (Wiener-Khinchin theorem) is the bridge between the time-domain view (autocorrelation) and the frequency-domain view (PSD)** — arguably the single most important fact in this whole document, since it's what lets you move freely between "how correlated is the process over time" and "how is its power distributed over frequency," which is exactly what's needed to reconcile the time-domain coherence-time picture and the frequency-domain coherence-bandwidth picture from the channel models notes.

### 4.2 Why "white" noise is called white

Taking the Fourier transform of the delta-function autocorrelation from Section 3.2:

```
S_N(f) = N0/2    (constant for all f)
```

**A flat PSD across all frequencies** — analogous to white light containing all colors (frequencies) with equal intensity, hence the name. This is the direct link between the time-domain definition (delta-function autocorrelation, Section 3.2) and the frequency-domain definition (flat PSD) — they are the *same* statement, just viewed through the Wiener-Khinchin transform. Interviewers sometimes ask "why is it called white noise" specifically to check you know this connection rather than just the flat-PSD fact in isolation.

### 4.3 PSD and the AWGN capacity formula — closing the loop with the capacity notes

The noise power `N = N0*B` used directly in the AWGN capacity formula (`C = B*log2(1+P/N)`, capacity notes Section 1.2) comes from **integrating the flat PSD `N0/2` over the (two-sided) bandwidth** occupied by the signal — this is literally where the `N0*B` noise-power term in the capacity formula comes from, connecting this section directly back to the very first document in this series.

### 4.4 PSD and filtering — a commonly tested relationship

If a WSS process `X(t)` with PSD `S_X(f)` passes through a linear time-invariant filter with frequency response `H(f)`, the output PSD is:

```
S_Y(f) = |H(f)|^2 * S_X(f)
```

**This is exactly the relationship underlying the matched filter's SNR-maximization result (detection theory notes, Section 2.3)** — the matched filter's optimality proof implicitly uses this PSD-shaping relationship to reason about how filtering reshapes the noise power spectrum relative to the signal.

### 4.5 PSD of fading processes — closing the loop with the channel models notes

The PSD of a time-varying fading process (as a function of Doppler frequency) is called the **Doppler power spectrum** (classically, the **Jakes spectrum** for isotropic scattering with a moving receiver). Its shape directly determines the coherence time via the inverse relationship between time-domain decorrelation and frequency-domain spread — this is the precise, formal version of "coherence time relates to Doppler spread" stated informally in the channel models notes (Section 4.4): coherence time and Doppler spread are a Fourier-transform pair, exactly analogous to how coherence bandwidth and delay spread (Section 4.5 of the channel models notes) are a Fourier-transform pair on the frequency/delay side. **Both of these "coherence X relates inversely to Y" facts from the channel models notes are really the same Wiener-Khinchin relationship, applied twice — once in time/Doppler, once in frequency/delay** — a strong, honest synthesis point if asked to connect the two.

---

## 5. Ergodicity (a commonly confused, adjacent concept — worth distinguishing clearly from stationarity)

### 5.1 Definition

A process is **ergodic** if time averages (computed from a single, sufficiently long realization) equal ensemble averages (computed by averaging across many independent realizations at a fixed time):

```
(1/T) * integral_0^T X(t) dt   ->   E[X(t)]    as T -> infinity   (mean-ergodic, for one example)
```

### 5.2 Why this distinction matters practically

In practice, you almost never have access to many independent realizations of a channel/noise process — you have **one** long recording. Ergodicity is the property that justifies estimating statistics (mean, autocorrelation, PSD) from time-averages of that single recording, and treating the result as if it were the true ensemble statistic. **Stationarity and ergodicity are logically distinct properties** — a process can be stationary but not ergodic (e.g., a process whose value is a single random constant per realization, drawn once at time zero and then simply held forever, is trivially stationary but clearly not ergodic — the time average just returns that one frozen realization's constant value, never converging to the true ensemble mean across different realizations). Worth having this or a similar minimal counterexample ready if asked to show you understand the two concepts aren't the same thing.

### 5.3 Connecting "ergodic" back to "ergodic capacity" — closing a loop from the very first document

This is worth stating explicitly, since it's an easy connection to miss under interview pressure: the term "**ergodic capacity**" (capacity notes, Section 2.2) is named for *exactly* this property — ergodic capacity is the rate achievable when the code is long enough that the transmitted codeword's performance, averaged over time, converges to the *ensemble*-averaged capacity `E_h[B*log2(1+gamma)]` — i.e., the coding scheme's time-average behavior over one long transmission matches the statistical (ensemble) average, which is precisely the definition of ergodicity given above. This naming connection is a nice, low-effort thing to mention if the conversation moves from probability theory into capacity — it shows the vocabulary isn't arbitrary.

---

## 6. Quick-reference formula sheet

```
WSS conditions:              mu_X(t) = mu_X (constant),  R_X(t1,t2) = R_X(t1-t2) = R_X(tau)
Autocorrelation at 0 lag:     R_X(0) = E[|X(t)|^2] = average power
White noise autocorrelation:  R_N(tau) = (N0/2)*delta(tau)
Wiener-Khinchin theorem:      S_X(f) = FourierTransform{ R_X(tau) }
White noise PSD:              S_N(f) = N0/2   (flat/constant)
Filtered process PSD:         S_Y(f) = |H(f)|^2 * S_X(f)
Ergodic mean:                 (1/T)*integral X(t)dt -> E[X(t)]  as T->infinity
```

---

## 7. Practice prompts

1. State the two conditions for wide-sense stationarity and explain why WSS (not SSS) is the standard working assumption in communications theory.
2. Explain why WSS implies SSS specifically for Gaussian processes, but not in general.
3. Explain, using the Wiener-Khinchin theorem, why "white" noise (delta-function autocorrelation) has a flat power spectral density — walk through the Fourier transform intuition, not just the result.
4. Derive/explain where the `N = N0*B` noise power term in the AWGN capacity formula comes from, starting from the flat noise PSD.
5. Explain the difference between stationarity and ergodicity, and give a simple example of a process that is stationary but not ergodic.
6. Explain why "ergodic capacity" is named the way it is — connect the probability-theory definition of ergodicity to the capacity-theory usage.
7. Explain how coherence time (time domain) and Doppler spread (frequency domain) are related via Wiener-Khinchin, and draw the parallel to how coherence bandwidth and delay spread are related the same way.
8. Connect this to your own research: snapshot-based DOA estimation methods (MUSIC, GS-SBL) rely on the array covariance matrix `R` being estimated from multiple snapshots — what stationarity/ergodicity assumption is implicitly being made when you replace the true covariance with a sample covariance averaged over snapshots?
