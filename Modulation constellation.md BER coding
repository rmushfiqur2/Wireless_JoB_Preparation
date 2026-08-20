# Modulation, Constellation Design, BER Curves, and Coding — Study Notes

Prepared as a companion to the Shannon capacity notes — this section covers how real systems approach the capacity limit: constellation design, error probability (BER) analysis, and channel coding.

---

## 1. Digital Modulation & Constellation Design

### 1.1 The basic idea

Shannon's theorem tells you the *rate limit*; modulation and coding tell you *how* to actually build a system that transmits bits over a physical channel. A **modulation scheme** maps groups of bits to points in a 2D **signal constellation** (in-phase `I` / quadrature `Q` plane), each point representing a distinct waveform (amplitude/phase combination) sent over the channel.

```
bits -> symbol mapping -> complex symbol s = I + jQ -> transmitted waveform
```

### 1.2 Common constellations (know these cold)

- **BPSK (Binary PSK):** 2 points on the real axis (`+-1`), 1 bit/symbol. Most robust to noise, lowest spectral efficiency.
- **QPSK (Quadrature PSK):** 4 points at 45°/135°/225°/315°, 2 bits/symbol. Equivalent to two independent BPSK streams on I and Q.
- **M-QAM (Quadrature Amplitude Modulation):** grid of `M` points varying in both amplitude and phase (e.g., 16-QAM = 4x4 grid, 64-QAM, 256-QAM, 1024-QAM). `log2(M)` bits/symbol.
- **M-PSK:** `M` points equally spaced on a circle (constant amplitude, varying phase only) — used when amplitude distortion (e.g., nonlinear power amplifiers) is a concern, since constant-envelope signals are more robust to amplifier nonlinearity.

**Key design tradeoff:** higher-order constellations (64-QAM, 256-QAM) pack more bits/symbol → higher spectral efficiency → but points are closer together → more sensitive to noise → need higher SNR to achieve the same error rate. This is the central tradeoff you should be ready to explain: **spectral efficiency vs. robustness**, and it's exactly why systems like LTE/5G NR use **adaptive modulation and coding (AMC)** — switching to lower-order constellations (QPSK) when SNR is poor (cell edge) and higher-order (256-QAM) when SNR is good (near the base station).

### 1.3 Constellation design principles

- **Maximize minimum Euclidean distance** `d_min` between constellation points for a given average energy constraint — this is the single design objective that drives error performance (see BER section below).
- **Gray coding:** assign bit patterns to adjacent constellation points so that neighboring points differ by only 1 bit. This matters because the most likely error (noise pushing a symbol to an adjacent point) then causes only a single bit error rather than multiple — dramatically improves the effective **bit** error rate even though the **symbol** error rate is unchanged.
- **Average energy normalization:** constellations are typically normalized so average symbol energy is fixed (e.g., `E_avg = 1`), so comparisons across schemes are fair — this is why higher-order QAM constellations look "denser" rather than simply "bigger."

---

## 2. Bit Error Rate (BER) Analysis

### 2.1 The basic detection problem

At the receiver, after passing through an AWGN channel, you observe `y = s + n` and must decide which constellation point `s` was sent. The optimal decision rule (for equally likely symbols, AWGN) is **minimum-distance detection** — pick the constellation point closest to `y` in Euclidean distance. This is the ML/MAP detector for this problem (ties directly into detection theory).

### 2.2 BER for BPSK (the baseline you should be able to derive/state [Detailed Math](https://claude.ai/chat/8d798066-bdcb-490c-a85c-6fbca03bc3f0))

```
P_b(BPSK) = Q( sqrt(2*Eb/N0) )
```

where `Q(x)` is the Gaussian tail function `Q(x) = (1/sqrt(2*pi)) * integral from x to infinity of exp(-t^2/2) dt`, and `Eb/N0` is energy-per-bit to noise-power-spectral-density ratio (the standard currency for comparing modulation/coding schemes independent of bandwidth/rate).

**Where this comes from (sketch):** the two BPSK symbols are separated by distance `d = 2*sqrt(Eb)`. An error occurs when noise pushes the received point past the midpoint (decision boundary) — i.e., when the noise component along the signal direction exceeds `sqrt(Eb)`. Since the noise is Gaussian with variance `N0/2`, this probability is exactly `Q(sqrt(2*Eb/N0))`.

### 2.3 BER for other schemes (know the trends, not necessarily every formula)

- **QPSK:** same BER as BPSK per bit — `P_b(QPSK) = Q(sqrt(2*Eb/N0))` — because QPSK is just two orthogonal BPSK streams (I and Q independently). QPSK gets you 2x the bits/symbol "for free" relative to BPSK bandwidth-wise, without a BER penalty, *when compared at the same Eb/N0*. This is a classic interview point: "why is QPSK often preferred over BPSK?" — same energy efficiency, double the rate.
- **M-QAM (approximate, high-SNR, Gray-coded):**
  ```
  P_b(M-QAM) ~ (4/log2(M)) * (1 - 1/sqrt(M)) * Q( sqrt( 3*log2(M)/(M-1) * Eb/N0 ) )
  ```
  You don't need to memorize this exactly, but you should understand the *shape*: as `M` grows, you need higher `Eb/N0` for the same BER — the argument of the `Q(.)` function shrinks as `M` grows (roughly as `1/M`), meaning the effective distance between points shrinks, so error probability rises unless you compensate with more power.
- **M-PSK:** similar flavor — BER degrades as `M` increases because points are packed more tightly around the circle (smaller angular separation → smaller Euclidean distance).

### 2.4 The BER curve — what it looks like and how to read it

A BER curve plots `P_b` (log scale, y-axis) vs. `Eb/N0` in dB (x-axis). Key features to recognize instantly:

- **Steep downward slope (waterfall region):** BER drops sharply once `Eb/N0` crosses a threshold — this threshold is the operating point systems are designed around.
- **Curves shift right as constellation order increases:** 64-QAM's curve sits to the right of QPSK's — same BER requires several dB more `Eb/N0`. Each doubling of bits/symbol via higher-order QAM typically costs roughly 3-4 dB in required SNR for the same target BER (a good ballpark number to quote).
- **Coding gain (see Section 3) shifts curves left:** a well-designed code lets you hit the same BER at *lower* `Eb/N0` — this leftward shift, measured in dB at a fixed BER, is literally the definition of "coding gain."
- **Diversity order shows up as slope:** in fading channels, the *slope* of the BER curve (on a log-log plot of BER vs. SNR) reflects the diversity order — higher diversity (more independent fading paths combined) gives a steeper slope, meaning BER falls off faster with increasing SNR. A single-antenna Rayleigh-fading link has BER falling only as `1/SNR` (much shallower than AWGN's exponential-like `Q()` falloff) — this is *why* fading is so damaging and why diversity techniques matter so much in practice.

---

## 3. Channel Coding

### 3.1 Why coding exists

Shannon's theorem says capacity `C = B*log2(1+SNR)` is achievable *with error probability approaching zero*, but the proof is non-constructive — it doesn't tell you *how*. Modulation alone (uncoded transmission) leaves you far from this limit. **Channel coding** adds structured redundancy so the receiver can detect and correct errors, pushing achievable performance toward the Shannon limit.

### 3.2 Coding gain — the key metric

**Coding gain** = the reduction in required `Eb/N0` (in dB) to achieve a target BER, compared to uncoded transmission. E.g., "this code provides 3 dB of coding gain at BER = 1e-5" means you need 3 dB less power to hit the same error rate.

- Comes at a cost: **redundancy** (code rate `r = k/n < 1`, sending `n` coded bits for every `k` information bits) means you either need more bandwidth for the same information rate, or you sacrifice some of your information rate for the same bandwidth — coding gain isn't free, it's a power-vs-bandwidth/rate tradeoff, conceptually similar to the modulation-order tradeoff above.

### 3.3 Families of codes (know the landscape, not implementation detail)

- **Block codes (e.g., Hamming, Reed-Solomon, BCH):** fixed-length codewords; classical algebraic constructions. Reed-Solomon still used for burst-error correction (e.g., DVB, storage).
- **Convolutional codes:** encode a continuous stream via a shift-register-based encoder; decoded with the **Viterbi algorithm** (maximum-likelihood sequence detection over the code trellis). Foundational in 2G/3G and still used as a building block/comparison point.
- **Turbo codes:** parallel concatenation of convolutional codes with an interleaver, decoded iteratively (two decoders "turbo" exchanging soft information). Used in 3G/4G LTE. Famous for getting within about 1 dB of the Shannon limit.
- **LDPC (Low-Density Parity-Check) codes:** sparse parity-check matrix, decoded with iterative **belief propagation / sum-product** message passing. Used for data channels in 5G NR (and Wi-Fi, DVB-S2). Near-capacity performance, highly parallelizable decoding (good for hardware).
- **Polar codes:** based on channel polarization (Arikan, 2008) — the only *provably* capacity-achieving code with an explicit, low-complexity construction as block length grows. Used for control channels in 5G NR (chosen for excellent short-blocklength performance, important for control signaling).

**Interview-relevant summary table:**

| Code | Decoding | Used in | Strength |
|---|---|---|---|
| Convolutional | Viterbi (ML sequence) | 2G/3G, legacy | Simple, well-understood |
| Turbo | Iterative (BCJR/MAP) | 3G/4G LTE data | ~1 dB from capacity |
| LDPC | Belief propagation | 5G NR data channel, Wi-Fi | Near-capacity, parallel decode |
| Polar | Successive cancellation (SC/SCL) | 5G NR control channel | Provably capacity-achieving, great at short blocklengths |

### 3.4 Soft vs. hard decision decoding

- **Hard-decision decoding:** receiver first makes a bit decision (0/1) for each received symbol, then decodes based on those hard bits. Simpler, but throws away information.
- **Soft-decision decoding:** receiver passes the raw (or likelihood-based) received values into the decoder, preserving confidence information. Soft decoding typically gains **about 2 dB** over hard-decision decoding for the same code — a number worth knowing, since it's a very standard interview fact.

### 3.5 Coded modulation: bringing it together

Real systems jointly design modulation and coding — hence "MCS" (**Modulation and Coding Scheme**) tables in LTE/NR, which pair a modulation order (QPSK/16-QAM/64-QAM/256-QAM) with a code rate, and adaptively select the best MCS index based on measured channel quality (CQI feedback) — directly implementing the AMC concept from Section 1.2. This is the practical, systems-level embodiment of "operate as close to the Shannon capacity curve as current SNR conditions allow."

---

## 4. Quick-reference formula sheet

```
BPSK BER:                 P_b = Q( sqrt(2*Eb/N0) )
QPSK BER:                 P_b = Q( sqrt(2*Eb/N0) )      [same as BPSK per bit]
M-QAM BER (approx):       P_b ~ (4/log2 M)(1 - 1/sqrt M) * Q( sqrt(3*log2(M)/(M-1) * Eb/N0) )
Code rate:                r = k/n   (k info bits, n coded bits)
Coding gain:               dB reduction in required Eb/N0 at fixed target BER, vs. uncoded
Soft vs hard decoding:     soft decoding gains ~2 dB over hard decoding
```

---

## 5. Practice prompts

1. Derive the BPSK BER formula from the minimum-distance decision rule and the Gaussian noise model.
2. Explain why QPSK achieves the same BER as BPSK at the same `Eb/N0`, despite carrying twice the bits/symbol.
3. Why does Gray coding matter for BER even though it doesn't change the symbol error rate?
4. Sketch (qualitatively) how a BER curve shifts as you move from QPSK -> 16-QAM -> 64-QAM at fixed transmit power, and explain the shift in terms of constellation geometry.
5. What is "coding gain," and how would you describe the fundamental tradeoff coding makes to achieve it?
6. Why does 5G NR use polar codes for control channels but LDPC for data channels? (Hint: think about blocklength regimes and reliability requirements.)
7. Connect this to your own research: sparse recovery / compressed sensing methods (like your GS-SBL work) also involve a form of "coding gain" intuition — redundancy/structure enabling recovery from noisy/incomplete measurements. How would you draw that analogy if asked?
