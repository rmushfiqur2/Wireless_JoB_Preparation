# Channel Models: Path Loss, Shadowing, Fading, and mmWave Blockage — Study Notes

Companion to the capacity and modulation/coding notes — this section covers how the wireless channel itself is modeled, from large-scale power decay down to small-scale fading and mmWave-specific effects.

---

## 1. The Three-Layer View of Channel Modeling

Wireless channel gain is conventionally decomposed into three multiplicative (additive in dB) components, separated by the *distance/time scale* over which they vary:

```
Received power (dB) = Transmit power - Path Loss - Shadowing + Small-scale fading (dB)
```

| Layer | Physical cause | Varies over |
|---|---|---|
| **Path loss** | Distance-dependent spreading/absorption of energy | Large scale (hundreds of meters / km) |
| **Shadowing** | Blocking by large obstacles (buildings, terrain) | Medium scale (tens of meters) |
| **Small-scale fading** | Multipath interference (constructive/destructive) | Small scale (fractions of a wavelength) |

This layered decomposition itself is a very common interview framing device — being able to draw this picture and place Rayleigh/Rician/path-loss/shadowing into it correctly is often worth more than any single formula.

---

## 2. Path Loss

### 2.1 Free-space path loss (Friis equation)

For a line-of-sight link in free space:

```
PL(dB) = 20*log10(d) + 20*log10(f) + 20*log10(4*pi/c)
```

or equivalently, received power via the **Friis transmission equation**:

```
Pr = Pt * Gt * Gr * (lambda / (4*pi*d))^2
```

`Pt`/`Pr` = transmit/receive power, `Gt`/`Gr` = antenna gains, `lambda` = wavelength, `d` = distance.

**Key scaling facts to have ready:**
- Path loss grows with the **square of distance** (`d^2` in linear power, i.e., `20*log10(d)` — a 6 dB loss per doubling of distance) in free space.
- Path loss grows with the **square of frequency** — higher frequency (e.g., mmWave) suffers more free-space loss for the same distance, all else equal. This is a big reason mmWave systems need large antenna arrays/beamforming to close the link budget.

### 2.2 Log-distance path loss model (more general, used for real environments)

```
PL(d)[dB] = PL(d0)[dB] + 10*n*log10(d/d0)
```

`d0` = reference distance, `n` = **path loss exponent**. `n=2` recovers free-space; realistic environments have `n` ranging roughly 2 (free space/open outdoor) to 4-6 (dense urban, indoor with obstructions, non-line-of-sight). The path loss exponent `n` is one of the most commonly cited "know this number" facts in wireless interviews.

### 2.3 LOS vs. NLOS

Standardized models (3GPP, ITU) typically define **separate path loss formulas for LOS and NLOS**, since NLOS paths (diffraction/reflection around obstacles) suffer significantly more loss than direct LOS paths at the same distance. Real systems (and simulators) often use a **LOS probability model** — a function of distance (and sometimes environment type) giving the probability a link is LOS vs. NLOS, then apply the corresponding path loss formula stochastically. This LOS-probability + separate-PL-curve approach is exactly the structure used in 3GPP's mmWave/UMa/UMi channel models.

---

## 3. Shadowing (Large-Scale/Slow Fading)

### 3.1 The model

Shadowing captures signal attenuation from large obstacles (buildings, hills, foliage) that block or partially block the LOS path. It's modeled as a **log-normal random variable** — i.e., Gaussian in the dB domain:

```
Received power (dB) = Pt - PL(d)[dB] - X_sigma
```

where `X_sigma ~ N(0, sigma^2)` is the shadowing term, with `sigma` typically in the range of **4-10 dB** depending on environment (lower for open/rural, higher for dense urban/indoor).

**Why log-normal, physically:** shadowing loss is the cumulative product of many independent multiplicative attenuation effects (each obstruction attenuates by some random factor); by the **Central Limit Theorem**, the *sum* of the log of these factors (i.e., the total loss in dB) tends toward Gaussian, hence "log-normal" in linear power.

### 3.2 Correlation

Shadowing is spatially correlated — two nearby locations experience similar shadowing (same buildings blocking the path), often modeled with an exponential spatial autocorrelation:

```
Corr(X_sigma(d1), X_sigma(d2)) = exp(-|d1-d2| / d_corr)
```

`d_corr` (correlation distance) is typically tens of meters. This matters for network planning/simulation (e.g., simulating handover, interference correlation between nearby users) — worth mentioning if the conversation goes toward system-level simulation.

---

## 4. Small-Scale Fading

### 4.1 Physical origin

Small-scale fading arises from **multipath propagation**: the receiver sees a sum of many delayed, phase-shifted copies of the transmitted signal (reflections, diffractions, scattering). These copies can add constructively or destructively depending on relative phase, causing rapid fluctuations in received signal amplitude over distances on the order of a wavelength.

### 4.2 Rayleigh fading

**When it applies:** many scattered paths, no dominant line-of-sight component (rich scattering, NLOS).

**Model:** by the Central Limit Theorem, with many independent multipath components, the received complex baseband signal's I and Q components are each approximately zero-mean Gaussian. The **envelope** `|h|` of a complex Gaussian is **Rayleigh-distributed**:

```
f(r) = (r/sigma^2) * exp(-r^2 / (2*sigma^2)),   r >= 0
```

and `|h|^2` (the instantaneous power/SNR-proportional quantity) is **exponentially distributed** — this exponential-power fact is exactly what underlies the "zero-outage capacity is zero" result from the capacity notes, since the exponential distribution has support all the way down to zero with positive density.

### 4.3 Rician fading

**When it applies:** a dominant LOS (or otherwise strong specular) path plus weaker scattered multipath.

**Model:** envelope follows a **Rician distribution**, parameterized by the **K-factor**:

```
K = (power in LOS component) / (power in scattered components)
```

- `K -> 0`: no LOS, Rician reduces to Rayleigh.
- `K -> infinity`: pure LOS, no fading (approaches AWGN — deterministic channel gain).
- Typical LOS mmWave/small-cell links: moderate-to-large `K` (e.g., 5-15 dB is common), reflecting a strong direct path plus some scattering.

**Interview framing:** be ready to state that Rician is the *general* model and Rayleigh/AWGN are its two limiting cases — a clean way to show you understand the whole family rather than three disconnected models.

### 4.4 Doppler spread / coherence time (the "how fast does it fade" axis)

Small-scale fading also varies **over time** due to mobility. Key quantities:

- **Doppler shift:** `f_d = v*cos(theta)/lambda`, max Doppler `f_d,max = v/lambda`.
- **Coherence time:** roughly `T_c ~ 1/f_d,max` — the time over which the channel stays approximately constant. This is exactly the quantity that determines whether you're in the "fast fading" (ergodic capacity, code spans many coherence times) or "slow/block fading" (outage capacity, code shorter than coherence time) regime from the capacity notes — worth explicitly connecting these two sections if asked.

### 4.5 Delay spread / coherence bandwidth (the frequency-domain analog)

- **Delay spread** `tau_rms`: spread of arrival times of multipath components.
- **Coherence bandwidth:** roughly `B_c ~ 1/tau_rms` — the bandwidth over which the channel's frequency response is approximately flat.
- **Flat fading vs. frequency-selective fading:** if signal bandwidth `<< B_c`, the channel is **flat fading** (all frequencies fade together). If signal bandwidth `>> B_c`, the channel is **frequency-selective** (different frequencies fade independently) — this is precisely why OFDM is attractive: it splits a wideband signal into many narrowband subcarriers, each experiencing approximately flat fading, turning a hard frequency-selective equalization problem into simple per-subcarrier scalar equalization.

---

## 5. mmWave-Specific Channel Behavior

mmWave (roughly 24-100+ GHz) channels differ enough from sub-6 GHz that they get their own modeling treatment. Key differences to know:

### 5.1 Higher path loss, compensated by large arrays

As shown in Section 2.1, free-space path loss scales with `f^2` — mmWave suffers substantially more path loss than sub-6 GHz at the same distance. This is compensated by **very large antenna arrays** (small wavelength → many antennas fit in a small area) providing high beamforming gain, which is why mmWave systems are inherently tied to massive MIMO/beamforming — you cannot separate the two topics in a real mmWave discussion.

### 5.2 Blockage — the defining mmWave challenge

mmWave signals are highly susceptible to **blockage**: they diffract poorly around obstacles (short wavelength → behaves more like light — sharp shadows) and are strongly attenuated by common materials (human body, foliage, even heavy rain at higher mmWave bands). This gives mmWave channels a qualitatively different character from sub-6 GHz:

- **Discrete LOS/NLOS/blocked states** rather than smooth fading — a mmWave link often behaves almost like a binary/tri-state process: strong LOS, weaker NLOS (via reflection), or effectively **outage-level blockage** (link nearly dead), rather than a continuously-varying fading envelope.
- **Human body blockage** is a specific, well-studied effect (e.g., a person walking between a UE and a small cell) — commonly modeled with a simple **knife-edge diffraction** model or empirically fit "blockage loss" of 15-40 dB, occurring on the timescale of pedestrian movement (much slower than small-scale multipath fading, much faster than large-scale shadowing — its own intermediate timescale).
- **Blockage probability models:** stochastic geometry approaches often model blockers (buildings, people) as a random spatial process (e.g., Boolean/Poisson line process) and derive the probability a link is blocked as a function of distance — this connects mmWave blockage modeling to the LOS-probability framework in Section 2.3, just with a more aggressive/sharper transition and additional blockage states.

### 5.3 Sparse multipath / clustered channel model

Unlike rich-scattering sub-6 GHz environments (which motivate the Rayleigh/CLT argument), mmWave channels typically have **few dominant multipath clusters** (the high path loss kills weak scattered paths before they reach the receiver with meaningful energy). This is commonly modeled with a **geometric/clustered channel model** (e.g., extended Saleh-Valenzuela model):

```
H = sum_{l=1}^{L} alpha_l * a_r(theta_r,l) * a_t(theta_t,l)^H
```

`L` = small number of paths/clusters, `alpha_l` = complex path gain, `a_r`/`a_t` = receive/transmit array steering vectors at angles `theta_r,l`/`theta_t,l`. This sparse, angularly-structured model is exactly the setting where **DOA estimation and sparse recovery methods (your GS-SBL background) become directly relevant** — mmWave channel estimation is often framed as a sparse recovery problem in the angular domain, because `L` is small relative to the array size. This is a strong, genuine bridge to your own research if the conversation goes there.

### 5.4 Frequent blockages -> reliability, not just fading, becomes the design driver

Because blockage events can be near-total outages rather than graceful fading dips, mmWave system design leans heavily on **multi-connectivity / beam-switching / relaying** strategies (rapidly switching beams or base stations when a blockage is detected) rather than relying purely on coding/diversity gain as in sub-6 GHz fading mitigation — worth mentioning as the systems-level consequence of the physical blockage model.

---

## 6. Quick-reference formula sheet

```
Friis (free space):        Pr = Pt*Gt*Gr*(lambda/(4*pi*d))^2
Free-space path loss:      PL(dB) = 20log10(d) + 20log10(f) + 20log10(4pi/c)
Log-distance path loss:    PL(d) = PL(d0) + 10*n*log10(d/d0)     [n ~ 2 (LOS) to 4-6 (NLOS/indoor)]
Shadowing:                 X_sigma ~ N(0, sigma^2) dB,  sigma ~ 4-10 dB
Rayleigh envelope:         f(r) = (r/sigma^2)*exp(-r^2/(2*sigma^2))
Rician K-factor:           K = P_LOS / P_scattered   (K->0: Rayleigh, K->inf: AWGN)
Coherence time:             T_c ~ 1/f_d,max
Coherence bandwidth:        B_c ~ 1/tau_rms
mmWave clustered channel:  H = sum_l alpha_l * a_r(theta_r,l) * a_t(theta_t,l)^H   [L small, sparse]
```

---

## 7. Practice prompts

1. Draw (mentally or on paper) the three-layer decomposition of received power and give the typical distance scale and physical cause of each layer.
2. Why does mmWave suffer more path loss than sub-6 GHz for the same distance, and what design choice directly compensates for it?
3. Explain why shadowing is modeled as log-normal, using a Central-Limit-Theorem argument.
4. State the relationship between Rayleigh and Rician fading in terms of the K-factor, and explain each limiting case physically.
5. Explain the connection between coherence time (channel modeling) and the ergodic-vs-outage capacity distinction (from the capacity notes) — this is a great "connect two sections" question a panel might ask.
6. Why is mmWave blockage often modeled as a near-binary LOS/NLOS/blocked state rather than smooth fading, and what system-level design responses does this motivate?
7. Connect this to your own research: explain why mmWave channel estimation is naturally framed as a sparse recovery problem, and how that relates to your GS-SBL / DOA estimation work.
