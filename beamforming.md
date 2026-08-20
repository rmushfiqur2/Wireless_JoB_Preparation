# Array Processing & Beamforming: DOA Estimation, MUSIC/ESPRIT, Beamforming Gain, MIMO Precoding — Study Notes

Companion to the estimation theory notes — this section is closest to your own dissertation work, so the framing here leans on connecting textbook material explicitly to your GS-SBL/DOA background rather than introducing it as new territory.

---

## 1. The Array Signal Model

### 1.1 Setup

Consider a uniform linear array (ULA) of `M` sensors/antennas receiving signals from `L` far-field sources at distinct angles `theta_1,...,theta_L`. The narrowband array output is:

```
y[n] = A(theta) * s[n] + w[n]
```

- `y[n]`: `M x 1` array snapshot at time `n`
- `A(theta) = [a(theta_1), ..., a(theta_L)]`: `M x L` **array manifold matrix**, columns are **steering vectors**
- `s[n]`: `L x 1` source signal amplitudes at time `n`
- `w[n]`: `M x 1` additive noise, typically modeled `~ N(0, sigma^2*I)`

### 1.2 The steering vector — the core object

For a ULA with element spacing `d`, wavelength `lambda`, the steering vector for angle `theta` is:

```
a(theta) = [1, exp(-j*2*pi*d/lambda*sin(theta)), ..., exp(-j*2*pi*(M-1)*d/lambda*sin(theta))]^T
```

**Physical meaning:** a plane wave arriving at angle `theta` hits each successive antenna element with a small extra propagation delay, producing a linear phase progression across the array. The steering vector encodes exactly this phase progression — it's the "fingerprint" of a source at angle `theta`. All of array processing (beamforming, DOA estimation) is fundamentally about exploiting this predictable phase structure.

**Half-wavelength spacing rule (`d = lambda/2`):** the standard choice to avoid spatial aliasing (analogous to the Nyquist sampling criterion, but in the spatial/angular domain rather than time/frequency). Worth having ready if asked "why lambda/2 spacing" — same math, different domain.

### 1.3 The covariance matrix — where DOA algorithms actually operate

Assuming sources and noise are uncorrelated, the array covariance matrix is:

```
R = E[y*y^H] = A * Rs * A^H + sigma^2 * I
```

`Rs = E[s*s^H]` is the `L x L` source covariance. In practice, `R` is estimated from `N` snapshots via the sample covariance:

```
R_hat = (1/N) * sum_{n=1}^{N} y[n] * y[n]^H
```

Essentially every classical DOA algorithm (beamforming, MUSIC, ESPRIT) works by extracting angle information from the eigenstructure of `R` (or `R_hat`) rather than from raw time-domain snapshots directly.

---

## 2. Conventional (Delay-and-Sum) Beamforming

### 2.1 The idea

The simplest DOA/beamforming approach: scan a candidate angle `theta`, apply the corresponding steering vector as weights, and measure output power:

```
P(theta) = a(theta)^H * R_hat * a(theta)
```

Peaks in `P(theta)` over a scan of candidate angles indicate likely source directions. Physically, this "steers" the array to add constructively for signals from `theta` (the array's phase-compensated outputs sum coherently) while signals from other angles add less coherently.

### 2.2 Limitations — why more sophisticated methods exist

- **Resolution limited by array aperture:** conventional beamforming's angular resolution is limited by the Rayleigh/beamwidth criterion, roughly `~lambda/(M*d)` — you cannot resolve two closely-spaced sources closer than this beamwidth, no matter the SNR. This is the classical **Rayleigh resolution limit**, and it's the key limitation that motivates every "high-resolution"/**superresolution** method below.
- Doesn't exploit the full statistical structure of the problem (doesn't distinguish signal-subspace vs. noise-subspace information) — leaves resolution on the table relative to what's theoretically achievable (recall the CRB from the estimation notes — conventional beamforming does not approach it in the closely-spaced-source regime).

---

## 3. MUSIC (MUltiple SIgnal Classification)

### 3.1 The key idea: subspace decomposition

MUSIC exploits an eigendecomposition of `R`:

```
R = U_s * Lambda_s * U_s^H + U_n * Lambda_n * U_n^H
```

- `U_s` (`M x L`): **signal subspace** eigenvectors, corresponding to the `L` largest eigenvalues.
- `U_n` (`M x (M-L)`): **noise subspace** eigenvectors, corresponding to the remaining `M-L` eigenvalues (theoretically all equal to `sigma^2` in the noiseless-model limit, since only `L` directions carry signal energy).

**The fundamental orthogonality fact MUSIC exploits:** the true steering vectors `a(theta_l)` for `l=1,...,L` lie *exactly* in the signal subspace, and are therefore **orthogonal to the entire noise subspace**:

```
a(theta_l)^H * U_n = 0   for the true source angles
```

### 3.2 The MUSIC pseudospectrum

Scan candidate angles and measure how close the steering vector is to orthogonal to the noise subspace:

```
P_MUSIC(theta) = 1 / ( a(theta)^H * U_n * U_n^H * a(theta) )
```

At the true source angles, the denominator goes to (near) zero, producing sharp peaks — in principle, **infinitely sharp** in the noiseless case, which is why MUSIC is called a **superresolution** method: unlike conventional beamforming, its resolution is not fundamentally limited by aperture the same way — with enough SNR/snapshots, it can resolve sources closer than the classical beamwidth limit.

### 3.3 MUSIC algorithm summary (be ready to state these steps in order)

1. Collect `N` snapshots, form sample covariance `R_hat`.
2. Eigendecompose `R_hat`; identify `L` (number of sources — often estimated separately, e.g. via information-theoretic criteria like AIC/MDL, or assumed known).
3. Partition eigenvectors into signal subspace `U_s` and noise subspace `U_n`.
4. Compute the MUSIC pseudospectrum `P_MUSIC(theta)` over a grid of candidate angles.
5. Source DOA estimates = locations of the `L` largest peaks.

### 3.4 Strengths and limitations (a favorite interview probe: "when would you NOT use MUSIC?")

- **Strength:** superresolution, works for arbitrary array geometries (not just ULA), well-understood asymptotic performance.
- **Limitation — requires knowing/estimating `L`:** performance degrades badly with an incorrect source-count estimate.
- **Limitation — needs enough snapshots for a good covariance estimate**, and **degrades sharply with correlated/coherent sources** (e.g., multipath scenarios where two "sources" are actually correlated reflections of the same signal) — the eigenvalue split between signal/noise subspace becomes ambiguous when sources are highly correlated, since correlated sources don't actually populate `L` distinct dimensions of the signal subspace as cleanly. (Spatial smoothing is a common fix, worth mentioning if this comes up.)
- **Computational cost:** requires a fine angular grid search — this grid-search cost, and MUSIC's sensitivity to model-order errors, is exactly the practical motivation for **sparse-recovery-based DOA methods (like your GS-SBL work)**, which reformulate DOA estimation as a sparse signal recovery problem over a discretized angular grid, naturally estimating the number of sources via sparsity rather than requiring it as a separate input, and tend to be more robust to correlated sources than classical subspace methods. This is a very natural, honest bridge to raise if MUSIC's limitations come up.

---

## 4. ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques)

### 4.1 The key idea: shift-invariance, no grid search

ESPRIT exploits a special structural property of ULAs (or more generally, any array with a **shift-invariant** structure — two identical, displaced subarrays): the array is split into two overlapping subarrays, and the phase relationship between corresponding signal-subspace components of the two subarrays directly encodes the DOA — **without any angular grid search**.

### 4.2 The rotational invariance property

For a ULA split into subarray 1 (elements 1 to M-1) and subarray 2 (elements 2 to M, i.e., shifted by one element), the steering vectors of the two subarrays are related by a simple diagonal phase matrix `Phi`:

```
A_2 = A_1 * Phi,     Phi = diag( exp(-j*2*pi*d/lambda*sin(theta_l)) )_{l=1..L}
```

The signal subspaces of the two subarrays are related by this same `Phi` (up to an invertible transform), and `Phi`'s eigenvalues can be recovered via a generalized eigenvalue problem on the signal-subspace matrices of the two subarrays — the eigenvalues directly give the phase shifts, hence the DOAs, in closed form:

```
theta_l = arcsin( -angle(Phi_l) * lambda / (2*pi*d) )
```

### 4.3 MUSIC vs. ESPRIT — the head-to-head comparison (very likely direct interview question)

| | MUSIC | ESPRIT |
|---|---|---|
| Requires grid search? | Yes (fine angular scan) | No (closed-form via eigenvalues) |
| Array geometry | Works for arbitrary (known) geometry | Requires shift-invariant structure (e.g., ULA, or paired subarrays) |
| Computational cost | Higher (grid search resolution vs. cost tradeoff) | Lower (one generalized eigenvalue problem) |
| Resolution | Superresolution, limited by grid fineness in implementation | Superresolution, exact (no grid quantization error) |
| Robustness to correlated sources | Degrades (needs spatial smoothing fix) | Also degrades, similar fixes apply |

**One-sentence summary to have ready:** *MUSIC trades computational cost for generality (any array geometry) via an explicit search, while ESPRIT trades generality for computational efficiency and exact (non-quantized) angle estimates by exploiting shift-invariant array structure.*

---

## 5. Beamforming Gain

### 5.1 The basic result

Combining `M` antenna outputs coherently (matched to a known/estimated steering vector, i.e., **maximum ratio combining** in the spatial domain) improves output SNR by a factor of `M` relative to a single antenna, assuming independent noise across elements:

```
SNR_array = M * SNR_single
```

**In dB:** beamforming gain `= 10*log10(M)` — e.g., doubling the number of antennas gives roughly 3 dB of gain. This is the array-processing analog of the matched-filter SNR result from the detection theory notes — same underlying idea (coherent combining maximizes SNR against white noise), just across spatial samples (antennas) instead of a time-domain waveform.

### 5.2 Beamforming vs. diversity vs. multiplexing — the conceptual map (ties back to the MIMO capacity notes)

- **Beamforming:** use array gain/directivity to boost SNR toward one direction (single stream) — trades array size for link-budget improvement, directly relevant to closing the mmWave link budget discussed in the channel models notes.
- **Spatial diversity:** use multiple antennas to combat fading (independent fading realizations combined for reliability) — relevant when the concern is robustness, not raw SNR.
- **Spatial multiplexing:** use multiple antennas to send **independent parallel data streams** (the MIMO capacity-scaling regime from the capacity notes) — trades array gain for higher raw data rate.

Be ready to state plainly: **you cannot simultaneously maximize all three** with a fixed number of antennas — this is the classical diversity-multiplexing-gain-tradeoff idea, worth restating here since beamforming is really "all the antenna gain devoted to a single, well-aimed direction," which is the opposite extreme from full spatial multiplexing.

---

## 6. MIMO Precoding — connecting array processing back to the MIMO capacity notes

### 6.1 The link to waterfilling

Recall from the capacity notes: with channel state at the transmitter (CSIT), MIMO capacity is maximized via **waterfilling power allocation across eigenmodes** of `H`. **Precoding** is the mechanism that actually implements this: the transmitter applies a precoding matrix `W` (built from the right singular vectors of `H`, via SVD `H = U*Sigma*V^H`) to project data streams onto the channel's eigenmodes before transmission, and the receiver applies a corresponding combining matrix (left singular vectors `U^H`) to separate the streams cleanly:

```
x = W * s   (precoding at transmitter, W built from V)
y_processed = U^H * y = Sigma * s + U^H*n   (post-combining: clean parallel channels)
```

This SVD-based precoding turns the full MIMO channel into a set of **independent, non-interfering scalar channels**, each with gain equal to a singular value of `H` — exactly the "parallel eigenmode" picture from the MIMO capacity notes, now shown as something you actually *build*, not just an abstraction for a capacity formula.

### 6.2 Codebook-based precoding (the practical LTE/NR version)

In practice (LTE/NR), full SVD-based precoding requires very accurate, high-overhead CSIT feedback. Real systems instead use **codebook-based precoding**: a finite set of predefined precoding matrices (indexed by a **Precoding Matrix Indicator, PMI**) is standardized; the receiver measures the channel, picks the codebook entry that best matches (e.g., maximizes achievable rate or minimizes distortion), and feeds back just the index — vastly reducing feedback overhead relative to sending the full channel matrix or precoder. Worth mentioning as the systems-level, standards-relevant version of the SVD precoding theory — a good way to show you can bridge textbook math to how Qualcomm's actual products/standards work operate.

### 6.3 Beamforming/precoding and DOA — the closing connection to your work

In mmWave/massive-MIMO systems, precoding is very often **angle-based** rather than full-CSIT SVD-based (channel estimation overhead is prohibitive with huge arrays) — the transmitter estimates dominant angles-of-departure (structurally the transmit-side analog of DOA/angle-of-*arrival* estimation) and beamforms toward them, exactly the **sparse, clustered, angular-domain channel model** from the mmWave channel notes. This is the cleanest possible closing link across this whole reading series: DOA/array processing (your core background) directly enables the angle-domain beamforming/precoding strategies that make mmWave/massive-MIMO systems practical — worth having as a one-paragraph "here's how my research connects to Qualcomm's systems" answer if asked to summarize your fit for the role.

---

## 7. Quick-reference formula sheet

```
Steering vector (ULA):      a(theta) = [1, e^{-j*2*pi*d/lambda*sin(theta)}, ..., e^{-j*2*pi*(M-1)*d/lambda*sin(theta)}]^T
Array covariance:            R = A*Rs*A^H + sigma^2*I
Conventional beamforming:    P(theta) = a(theta)^H * R_hat * a(theta)
MUSIC pseudospectrum:        P_MUSIC(theta) = 1 / (a(theta)^H * U_n * U_n^H * a(theta))
ESPRIT invariance:           A_2 = A_1 * Phi,  theta_l from angle(eig(Phi))
Beamforming gain:             SNR_array = M * SNR_single   (~ 10*log10(M) dB)
MIMO precoding (SVD):        H = U*Sigma*V^H,  x = V*s,  y_processed = U^H*y = Sigma*s + noise
```

---

## 8. Practice prompts

1. Derive/explain why the true source steering vectors are orthogonal to the noise subspace in the noiseless array model — this is the one-sentence core of why MUSIC works.
2. Walk through the MUSIC algorithm steps from raw snapshots to DOA estimates, out loud, without notes.
3. Explain ESPRIT's rotational invariance property and why it avoids a grid search — contrast directly with MUSIC.
4. State the Rayleigh resolution limit for conventional beamforming and explain, conceptually, why MUSIC/ESPRIT can beat it (superresolution).
5. Explain why MUSIC degrades for highly correlated/coherent sources, and name at least one practical fix.
6. Connect SVD-based MIMO precoding to the waterfilling result from the capacity notes — explain how they're the same underlying idea expressed two different ways (capacity allocation vs. actual signal processing).
7. Explain, in your own words, how DOA/angle estimation connects to mmWave beamforming and precoding — this is close to a direct summary of your dissertation's relevance to Qualcomm's systems work, so have a tight, confident version ready.
8. Connect this to your own research directly: how does your GS-SBL approach compare to MUSIC and ESPRIT in terms of assumptions (correlated sources, known model order, grid search), and what does it gain by framing DOA estimation as sparse Bayesian recovery instead?
