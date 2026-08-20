# Linear Algebra: Eigendecomposition, SVD, Matrix Rank — Study Notes

Companion to the array/beamforming and MIMO capacity notes — this section is the pure-math machinery underneath nearly everything covered so far (MUSIC's subspace split, MIMO waterfilling, compressed sensing). Framed here explicitly around *where each concept already showed up* in earlier sections, since that's exactly how a panel is likely to probe it — not in isolation, but as "why does this matrix have that property."

---

## 1. Eigendecomposition

### 1.1 Definition

For a square matrix `A` (n x n), an eigenvector `v` and eigenvalue `lambda` satisfy:

```
A*v = lambda*v
```

`v` is a direction that `A` only *stretches* (by factor `lambda`), never rotates. If `A` has `n` linearly independent eigenvectors, it can be **diagonalized**:

```
A = V * Lambda * V^{-1}
```

`V` = matrix of eigenvectors (columns), `Lambda` = diagonal matrix of eigenvalues.

### 1.2 The special, important case: symmetric/Hermitian matrices

**Every covariance matrix in this whole prep series (`R` in the MUSIC notes, `H*H^H` in the MIMO capacity notes) is Hermitian** (`A = A^H`, i.e., `A` equals its own conjugate transpose — the real-valued special case is "symmetric," `A = A^T`). Hermitian matrices have three properties that make them the backbone of everything above, and you should be able to state all three cold:

1. **Eigenvalues are real** (even if `A` itself is complex) — this is *why* it makes sense to talk about a covariance matrix's eigenvalues as "power levels" (Section 3.2, MIMO notes) — they're guaranteed to be real, non-negative numbers, not complex.
2. **Eigenvectors are orthogonal** (can be chosen orthonormal) — this is *exactly why* MUSIC's signal/noise subspace split (Section 3.1, array notes) is a clean orthogonal decomposition, not an approximate one.
3. **Diagonalization becomes `A = V*Lambda*V^H`** (i.e., `V^{-1} = V^H` since `V` is unitary/orthonormal) — no matrix inversion needed, just a conjugate transpose. This is called the **spectral theorem**.

### 1.3 Positive semi-definite (PSD) matrices

A Hermitian matrix `A` is **positive semi-definite** if `x^H*A*x >= 0` for all vectors `x`. Covariance matrices are always PSD (variance can't be negative) — this is *why* eigenvalues of a covariance matrix are guaranteed non-negative (property 1 above, sharpened): PSD-ness is exactly the condition that rules out negative eigenvalues.

### 1.4 Where eigendecomposition already appeared in this prep series

- **MIMO capacity (capacity notes, Section 3.2):** `H*H^H` is Hermitian PSD; its eigenvalues `lambda_i` are literally the per-eigenmode "channel gains" that appear inside each `log2(1+...)` term.
- **MUSIC (array notes, Section 3.1):** the covariance `R`'s eigendecomposition directly splits into signal subspace (large eigenvalues) vs. noise subspace (eigenvalues near `sigma^2`) — this split is only clean and orthogonal *because* `R` is Hermitian.
- **Waterfilling (capacity notes, Section 3.3):** the "buckets" being filled are literally the eigenvalues `lambda_i` of `H*H^H`.

---

## 2. Singular Value Decomposition (SVD)

### 2.1 Why you need SVD, not just eigendecomposition

Eigendecomposition only applies to **square** matrices, and only diagonalizes nicely when eigenvectors are guaranteed orthogonal (Hermitian case). But `H` (the MIMO channel matrix) is generally **rectangular** (`Nr x Nt`, and `Nr != Nt` in general) and not itself Hermitian. **SVD is the generalization that works for any matrix, square or not.**

### 2.2 Definition

Any matrix `H` (`m x n`) can be decomposed as:

```
H = U * Sigma * V^H
```

- `U` (`m x m`): orthonormal (**left singular vectors**)
- `V` (`n x n`): orthonormal (**right singular vectors**)
- `Sigma` (`m x n`): diagonal (in the leading `min(m,n) x min(m,n)` block) with **non-negative real singular values** `sigma_1 >= sigma_2 >= ... >= 0` on the diagonal

### 2.3 The direct link between SVD and eigendecomposition (a very likely direct interview question)

```
H^H*H = V * Sigma^H*Sigma * V^H     ->   eigenvalues of H^H*H are sigma_i^2, eigenvectors are V
H*H^H = U * Sigma*Sigma^H * U^H     ->   eigenvalues of H*H^H are sigma_i^2, eigenvectors are U
```

**One-sentence summary to have ready:** *the singular values of `H` are the square roots of the eigenvalues of `H^H*H` (or `H*H^H`) — SVD is essentially "eigendecomposition applied to the two Hermitian PSD matrices you can build from a rectangular `H`," and the two sets of singular vectors (`U`, `V`) are the eigenvectors of those two matrices respectively.*

### 2.4 Where SVD already appeared in this prep series

- **MIMO capacity notes (Section 3.2):** `lambda_i` (used directly in the capacity formula `sum log2(1+...)`) is literally `sigma_i^2` for `H` — the notes used `H*H^H`'s eigenvalues, which by the identity above are exactly the squared singular values of `H`.
- **MIMO precoding (array/beamforming notes, Section 6.1):** SVD-based precoding uses `V` (right singular vectors) as the precoding matrix and `U^H` as the receive combiner — this *is* the practical, "build a real system" use of SVD: it turns `H` into a clean diagonal `Sigma` matrix, i.e., independent parallel scalar channels.

### 2.5 Low-rank approximation (Eckart-Young theorem) — a useful, often-asked SVD fact

Truncating the SVD to the top `k` singular values/vectors gives the **best rank-`k` approximation** to `H` in the least-squares (Frobenius norm) sense:

```
H_k = sum_{i=1}^{k} sigma_i * u_i * v_i^H
```

This is the mathematical basis for PCA, compression, and (relevantly) **reduced-rank channel estimation/feedback** — if a channel matrix is effectively low-rank (few dominant eigenmodes/paths, as in the sparse mmWave clustered model from the channel-models notes), you only need to feedback/estimate the top few singular vectors, not the full matrix. Another clean cross-section link: mmWave's sparse angular structure (few dominant clusters) is *why* mmWave channel matrices are effectively low-rank, which is *why* SVD-based dimensionality reduction is so effective there.

---

## 3. Matrix Rank

### 3.1 Definition

`rank(A)` = the number of linearly independent rows (equivalently, columns) of `A` = the number of **non-zero singular values** of `A`. This last equivalence is the cleanest working definition to use once you have SVD in hand — it directly connects rank to something computable and already discussed above.

### 3.2 Rank and MIMO capacity — the direct link (ties Sections 2-3 of the capacity notes together)

Recall from the capacity notes: MIMO capacity scales with `min(Nt,Nr)` non-trivial eigenmodes. More precisely: **the number of usable parallel spatial channels is exactly `rank(H)`**, not `min(Nt,Nr)` in general — `min(Nt,Nr)` is only an *upper bound* on rank (a matrix can be rank-deficient even if it's "big enough" dimensionally). This is an important sharpening to have ready: a `4x4` MIMO channel does not automatically give you 4 parallel streams — it gives you `rank(H) <= 4` streams, and `rank(H)` depends on the actual propagation environment.

### 3.3 What makes `rank(H)` small in practice — ties directly to the channel models notes

- **Strong LOS / poor scattering:** if the environment has a single dominant path (e.g., strong LOS with little multipath), `H` is close to **rank-1** — effectively only one usable spatial stream regardless of how many antennas you have. This is exactly the point flagged in the MIMO capacity notes (Section 3.4) about why real-world MIMO gains depend on a scattering-rich environment.
- **Highly correlated antenna elements:** if antennas are spaced too closely (or in an environment with very limited angular spread), the rows/columns of `H` become nearly linearly dependent — same effect, different physical cause.
- **Correlated sources in the array-processing notes:** this is the exact same underlying phenomenon as MUSIC's "coherent sources" problem (array notes, Section 3.4) — coherent/correlated sources make the *source covariance* `Rs` rank-deficient, which is why the signal-subspace/noise-subspace eigenvalue split becomes ambiguous. **Rank deficiency is the single mathematical thread connecting "MUSIC fails on correlated sources" and "MIMO capacity gain requires scattering-rich channels"** — genuinely worth stating explicitly if you get the chance, since it demonstrates you see the same underlying math showing up in two "different" topics.

### 3.4 Rank and compressed sensing / sparse recovery — the bridge to your own work

Sparse recovery problems (like GS-SBL) are often phrased in terms of a measurement matrix's rank/null-space properties: recovering a sparse vector `x` from `y = Phi*x` (fewer measurements than unknowns, so `Phi` is "wide" and clearly rank-deficient/underdetermined in the ordinary linear-algebra sense) is only possible because of the *extra* assumption of sparsity — ordinary rank arguments alone would say the problem is unsolvable (infinitely many `x` satisfy `y=Phi*x` when `Phi` is rank-deficient relative to the number of unknowns), and sparsity is exactly the side information that resolves this ambiguity. This is a clean way to explain, from pure linear algebra, *why* compressed sensing needs to go beyond classical rank/least-squares arguments — a good one-paragraph answer if someone asks "why can't you just solve this with a pseudo-inverse?"

---

## 4. Quick-reference formula sheet

```
Eigendecomposition:           A*v = lambda*v,   A = V*Lambda*V^{-1}   (V*Lambda*V^H if A Hermitian)
Hermitian properties:          real eigenvalues, orthogonal eigenvectors, A = V*Lambda*V^H
PSD condition:                 x^H*A*x >= 0 for all x   <=>  all eigenvalues >= 0
SVD:                          H = U*Sigma*V^H   (any m x n matrix)
SVD <-> eigendecomposition:    sigma_i^2 = eigenvalues of H^H*H (and of H*H^H)
Low-rank approx (Eckart-Young): H_k = sum_{i=1}^k sigma_i*u_i*v_i^H  (best rank-k approx, Frobenius norm)
Rank:                         rank(A) = # linearly independent rows/cols = # non-zero singular values
```

---

## 5. Practice prompts

1. State the three key properties of Hermitian matrices (real eigenvalues, orthogonal eigenvectors, `V^{-1}=V^H`) and explain why each one matters for MUSIC's subspace decomposition specifically.
2. Derive the relationship between the singular values of `H` and the eigenvalues of `H^H*H` — show the algebra, don't just state it.
3. Explain why MIMO capacity depends on `rank(H)`, not just `min(Nt,Nr)`, and give a physical scenario (from the channel models notes) where `rank(H)` would be much less than `min(Nt,Nr)`.
4. Explain the Eckart-Young low-rank approximation result in one sentence, and connect it to why mmWave channel matrices (sparse angular structure) are good candidates for low-rank compression.
5. Explain, using rank/null-space language, why an underdetermined linear system `y=Phi*x` is generally unsolvable without extra structure — and why sparsity is the specific structure that rescues compressed sensing/GS-SBL-style recovery.
6. State the single mathematical concept (rank deficiency) that connects "MUSIC degrades with correlated sources" and "MIMO capacity gain needs scattering-rich channels" — practice saying this connection out loud in under 30 seconds, since it's exactly the kind of cross-topic synthesis a panel listens for.
