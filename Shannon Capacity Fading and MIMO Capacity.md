# Shannon Capacity, Fading, and MIMO Capacity — Study Notes

Prepared for interview prep on communications theory fundamentals: AWGN capacity, ergodic vs. outage capacity under fading, and MIMO capacity with waterfilling.

---

## 1. AWGN Channel Capacity (Shannon-Hartley)

### 1.1 The setup

Consider the simplest channel model:

```
y = x + n
```

where `x` is the transmitted signal (power constraint `E[x^2] <= P`), `n` is additive white Gaussian noise with power spectral density `N0/2`, and the channel has bandwidth `B`. The noise power in bandwidth `B` is `N = N0 * B`.

### 1.2 The capacity formula

Shannon's theorem says the capacity (maximum error-free data rate) is:

```
C = B * log2(1 + P/N)  bits/second
```

`P/N` is the **SNR**. This is the single most important formula in the whole subject — you should be able to write it blind and explain every symbol.

**Key intuition:** capacity grows *logarithmically* with SNR but *linearly* with bandwidth. This is why modern systems (5G, Wi-Fi 6/7) chase bandwidth (wider channels, mmWave, carrier aggregation) rather than just cranking transmit power — doubling power only adds `log2(2) = 1 bit/s/Hz`, while doubling bandwidth doubles capacity outright.

### 1.3 Where the formula comes from (sketch derivation — be ready to reproduce this)

1. Model the channel as a discrete-time AWGN channel with `2B` real dimensions per second (Nyquist sampling of a band-limited signal).
2. Capacity of a real AWGN channel (single use) with input power constraint `P` and noise variance `N` is:
   ```
   C_per_use = (1/2) * log2(1 + P/N)   bits/channel use
   ```
   This comes from maximizing mutual information `I(X;Y) = h(Y) - h(Y|X)` over all input distributions with `E[X^2] <= P`. The Gaussian input maximizes differential entropy `h(Y)` for a fixed second moment, giving `h(Y) = (1/2)log2(2*pi*e*(P+N))` and `h(Y|X) = h(N) = (1/2)log2(2*pi*e*N)`. Subtracting gives the `(1/2)log2(1+P/N)` result. [Detailed Math](https://claude.ai/chat/ba933bf0-2a16-4a38-9028-f08e24a4e376)
3. Multiply by `2B` channel uses per second → the factor of 2 cancels the 1/2, giving `C = B*log2(1+P/N)`.

**Why the Gaussian input is optimal:** for a fixed variance, the Gaussian distribution maximizes differential entropy among all real-valued distributions. Since noise is already Gaussian, Gaussian signaling matches the "worst-case" noise statistics and maximizes the achievable rate. This is a classic result worth being able to state even if you don't reproduce the full max-entropy argument.

### 1.4 Things interviewers like to probe

- **Bandwidth-limited vs. power-limited regime:** at high SNR, `C ~ B*log2(SNR)` (power-limited gains are small, add bandwidth). At low SNR, `log2(1+x) ~ x/ln(2)`, so `C ~ B*P/(N0*B*ln2) = P/(N0*ln2)` — capacity becomes *independent of bandwidth* and linear in power. This is the regime relevant to deep-space/low-power links.
- **Spectral efficiency:** `C/B = log2(1+SNR)` in bits/s/Hz — this is the number practical systems (LTE, NR) are benchmarked against; real systems achieve some fraction of this (a "gap to capacity") due to practical coding/modulation constraints.
- **Capacity vs. rate achieved by practical schemes:** know that QAM/coded systems approach but never exceed the Shannon limit; LDPC/turbo/polar codes get within a fraction of a dB.

---

## 2. Fading Channels: Ergodic vs. Outage Capacity

Real wireless channels aren't static AWGN — the channel gain `h` varies randomly over time (multipath fading, shadowing). This changes what "capacity" even means, because now there are two different useful notions.

### 2.1 The fading channel model

```
y = h*x + n
```

`h` is the (possibly complex) channel gain, often modeled as Rayleigh-distributed in magnitude (no line-of-sight) or Rician (with LOS component). The **instantaneous SNR** is:

```
gamma = |h|^2 * P / N
```

which is now a random variable, since `h` is random.

### 2.2 Ergodic capacity

**Definition:** the long-run average rate achievable when the channel varies over many independent fading states, and the code spans enough time to "see" the full distribution of `h` (i.e., the code is long relative to the channel's coherence time).

```
C_ergodic = E_h[ B * log2(1 + gamma) ]
```

You average the instantaneous Shannon capacity over the fading distribution. This is the right notion of capacity for **delay-tolerant** traffic — the transmitter doesn't need every single fading realization to support the rate, only the *average* over time.

**Key property:** by Jensen's inequality (since `log2(1+x)` is concave), 
```
C_ergodic <= B*log2(1 + E[gamma])
```
Fading *reduces* average capacity relative to an AWGN channel with the same average SNR — variability in SNR hurts you more than the average helps you (concavity argument), unless you can adapt power/rate to the channel.

### 2.3 Outage capacity

**Definition:** relevant when the code is short relative to coherence time (channel is roughly *constant* over the codeword — "slow fading" / block fading). In this regime, you can't average over the fading — you're stuck with whatever `h` you get for that block. If the instantaneous channel can't support your target rate `R`, the transmission is in **outage**.

```
P_out(R) = Pr[ B*log2(1+gamma) < R ]
```

The **outage capacity** `C_out(epsilon)` is the largest rate `R` such that `P_out(R) <= epsilon`, for a target outage probability `epsilon` (e.g., 1%).

This is the right notion for **delay-constrained** traffic (voice, real-time control) where you can't wait many coherence blocks to average things out — you pick a rate, accept that some fraction `epsilon` of blocks will fail, and that's your operating point.

### 2.4 Ergodic vs. outage — the key contrast (a favorite interview question)

| | Ergodic capacity | Outage capacity |
|---|---|---|
| Applies when | Fast fading / long codes spanning many coherence blocks | Slow (block) fading / short codes, one coherence block |
| What you compute | Expectation of rate over fading distribution | Probability that instantaneous rate falls below target |
| Traffic type | Delay-tolerant | Delay-constrained |
| Behavior at low outage target | N/A | Rate → 0 as `epsilon → 0` in Rayleigh fading (zero-outage capacity is 0!) — because `P(gamma near 0)` is always nonzero for Rayleigh |

**The "zero-outage capacity is zero" fact** for Rayleigh fading is a classic gotcha: since Rayleigh fading has support all the way down to `|h|=0`, no positive rate can be guaranteed with zero probability of outage. This is exactly why practical systems tolerate some outage probability rather than demanding perfect reliability.

### 2.5 Diversity as the fix

Since fading hurts capacity/reliability, systems fight it with **diversity** — multiple independent looks at the channel (time, frequency, space/antennas) so the probability that *all* paths are in deep fade simultaneously drops sharply. This is the conceptual bridge into MIMO (see below) and into why array/DOA-style processing (your own background) is valuable: multiple antennas offer spatial diversity and/or the ability to combine signals coherently (beamforming gain).

### 2.6 Common fading distributions (know these cold)

- **Rayleigh:** `|h|` Rayleigh-distributed, `|h|^2` exponentially distributed. No LOS component. Worst-case fading in terms of the "zero-outage capacity is zero" property.
- **Rician:** LOS + scattered multipath; parameterized by `K`-factor (ratio of LOS power to scattered power). As `K -> infinity`, Rician → AWGN (no fading). As `K -> 0`, Rician → Rayleigh.
- **Log-normal shadowing:** large-scale variation (buildings, terrain) — captured with the channel gain's dB value being Gaussian-distributed. Often layered on top of Rayleigh/Rician small-scale fading.

---

## 3. MIMO Capacity

### 3.1 The MIMO channel model

With `Nt` transmit and `Nr` receive antennas:

```
y = H*x + n
```

`H` is the `Nr x Nt` channel matrix, `x` is the `Nt x 1` transmit vector, `n` is `Nr x 1` AWGN with covariance `N0*I`.

### 3.2 Capacity with channel known at receiver only (no CSIT)

If the transmitter doesn't know `H` (only the receiver does), the natural strategy is to split power equally across transmit antennas: `x`'s covariance is `(P/Nt)*I`. Capacity is [Detailed Math](https://claude.ai/chat/1e4ba849-deb7-429a-9ef8-8501fab20274):

```
C = log2( det( I + (P/(Nt*N0)) * H * H^H ) )
```

Using the eigenvalue decomposition of `H*H^H` (eigenvalues `lambda_i`, `i=1..r`, `r = rank(H)`):

```
C = sum_i log2( 1 + (P/Nt)*lambda_i / N0 )
```

**Intuition:** MIMO turns one noisy scalar channel into `r` parallel independent "eigenmode" sub-channels, each behaving like its own scalar AWGN channel with gain `lambda_i`. This is the origin of the famous **MIMO capacity scales linearly with min(Nt, Nr)** result (at high SNR, roughly `min(Nt,Nr)` non-trivial eigenmodes each contribute a `log2(SNR)` term).

### 3.3 Capacity with channel known at transmitter (CSIT): waterfilling

If the transmitter *does* know `H` (e.g., via feedback), it can do better than equal power allocation by optimizing how much power to put on each eigenmode. This is the **MIMO waterfilling** problem:

```
maximize   sum_i log2( 1 + p_i * lambda_i / N0 )
subject to sum_i p_i <= P,  p_i >= 0
```

Using Lagrangian/KKT conditions, the solution is the classic **waterfilling** allocation:

```
p_i = max( mu - N0/lambda_i , 0 )
```

where `mu` (the "water level") is chosen so `sum_i p_i = P`.

**Intuition (the "pouring water" picture):** think of `N0/lambda_i` as the "floor height" of bucket `i` — weaker eigenmodes (small `lambda_i`) have a higher floor. You pour a fixed volume of water (total power `P`) over all the buckets; it settles at a common level `mu`. Strong eigenmodes (low floor) get filled deep; weak eigenmodes below the waterline get filled shallow; eigenmodes with floor above `mu` get **zero power** — you simply don't use very weak channels at all.

This is the same waterfilling principle that appears in frequency-selective OFDM channels (allocate more power to strong subcarriers, none to very weak ones) — same math, different context (spatial eigenmodes vs. frequency subcarriers).

### 3.4 Why this matters for a systems interview

- **Diversity-multiplexing tradeoff:** MIMO can be used for either spatial multiplexing (parallel data streams → more capacity) or diversity/beamforming (more reliability) — you generally cannot maximize both simultaneously; be ready to explain this tradeoff conceptually.
- **Rank of H matters:** in strong LOS conditions, `H` can be near rank-1 (all antennas see basically the same path), which kills the multiplexing gain even with many antennas — this is why real-world MIMO gains depend heavily on scattering-rich environments.
- **Massive MIMO / Giga-MIMO (Qualcomm's current 6G push):** scaling `Nt`/`Nr` way up increases available eigenmodes and lets beamforming concentrate energy — worth connecting your DOA/array-processing background here, since beamforming and MIMO precoding are deeply related to direction-of-arrival estimation.

---

## 4. Quick-reference formula sheet

```
AWGN capacity:            C = B * log2(1 + P/N)
Spectral efficiency:      C/B = log2(1 + SNR)   [bits/s/Hz]
Ergodic (fading):         C_ergodic = E_h[ B*log2(1+gamma) ]
Outage probability:       P_out(R) = Pr[ B*log2(1+gamma) < R ]
MIMO capacity (no CSIT):  C = sum_i log2(1 + (P/Nt)*lambda_i/N0)
Waterfilling (CSIT):      p_i = max(mu - N0/lambda_i, 0), sum p_i = P
```

---

## 5. Practice prompts (try deriving/answering without notes)

1. Derive the AWGN capacity formula from mutual information, starting from `I(X;Y) = h(Y) - h(Y|X)`.
2. Explain, using Jensen's inequality, why `C_ergodic <= C_AWGN` at the same average SNR.
3. Why is the zero-outage capacity of a Rayleigh fading channel exactly zero? What changes for Rician fading with large `K`?
4. Sketch the waterfilling solution for a 2x2 MIMO channel with eigenvalues `lambda_1 = 4`, `lambda_2 = 1` (pick a numeric power budget and actually solve for `p_1, p_2`).
5. Explain in one or two sentences why MIMO capacity scaling is often summarized as "linear in `min(Nt,Nr)`" rather than in `Nt*Nr`.
6. Connect this to your own research: how does DOA estimation / array processing relate to the eigenstructure of `H` in MIMO systems?
