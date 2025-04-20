# Point process

# Part 1 – What *Is* a Point Process?

## 1. Intuition
A **point process** is a probabilistic mechanism that scatters dimension‑less events in space—or in space × time—such that every realization is a dot map.

## 2. Motivation
* Goes beyond aggregate counts: it captures **where** and **when** events occur.
* Separates first‑order inhomogeneity (trend) from second‑order interaction (clustering or inhibition).

## 3. Formal Definition
Let the study window be $W \subset \mathbb{R}^2$ and the observation interval be $T \subset \mathbb{R}$.  
A realization is  
$$\mathbf{X} = \{(s_i, t_i)\}_{i=1}^{N} \subset W \times T,$$
where $(s_i, t_i)$ denotes the spatial location and timestamp of the $i$‑th event.  
Randomness enters through both **how many** points appear and **where/when** they appear.

## 4. First‑Order Property
The **intensity function** is defined as
$$\lambda(s, t) = \lim_{\lvert \Delta s \rvert,\, \Delta t \to 0} 
\frac{\mathbb{E}\,[N(\Delta s, \Delta t)]}
     {\lvert \Delta s \rvert \, \Delta t},$$
representing the expected number of events per unit area and per unit time at $(s, t)$.

* **Homogeneous** process: $\lambda(s, t)$ is constant.  
* **Inhomogeneous** process: $\lambda(s, t)$ varies with space, time, or covariates.

# Part 2 – Complete Spatial Randomness & the K‑Function

## 1. Intuition
**Complete Spatial Randomness (CSR)** posits that events occur independently and uniformly across the study region.  
* If CSR holds, any apparent clusters are purely accidental.  
* Rejecting CSR is the first step toward more realistic models (e.g., inhomogeneous Poisson, Cox, Hawkes).

## 2. Ripley’s $K$‑Function
For a 2‑D stationary point process with intensity $\lambda$,
\[
K(r) \;=\; \frac{1}{\lambda}\,
\mathbb{E}\bigl[\text{number of extra points within distance } r
\text{ of an arbitrary point}\bigr].
\]

* Under CSR (homogeneous Poisson) in the plane,
  \[
  K_{\text{CSR}}(r) = \pi r^{2}.
  \]
* A **pair‑correlation function** $g(r) = \dfrac{\mathrm{d}K(r)}{2\pi r\,\mathrm{d}r}$ gives local interaction:
  * $g(r) > 1$ ⇒ clustering at scale $r$  
  * $g(r) < 1$ ⇒ inhibition / regularity at scale $r$

## 3. Edge Correction
Points near the study window boundary have their neighborhoods truncated.  
Common corrections:
| Method | Idea | Use‑case |
|--------|------|----------|
| Border | Discard reference points closer than $r$ to the edge | Fast, biased if window small |
| Ripley isotropic | Weight each pair by reciprocal of the proportion of circle inside window | Default in `spatstat` |
| Translation | Translate edges onto themselves | Good for rectangular windows |

## 4. Monte Carlo CSR Test (Global Envelope)
1. Compute the empirical $\hat{K}(r)$ from data with edge correction.  
2. Simulate $M$ CSR patterns with the same intensity and window.  
3. For each simulation $m$, calculate $\hat{K}^{(m)}(r)$.  
4. Build upper/lower envelopes, e.g. 2.5⁠–⁠97.5 percentiles.  
5. **Reject CSR** if the observed $\hat{K}(r)$ leaves the envelope for any $r$.

> **Interpretation tip**  
> Plot the **$L$‑function** $L(r)=\sqrt{\hat{K}(r)/\pi}$; CSR corresponds to $L(r)=r$, so deviations are easier to see.

## 5. Example (R / `spatstat`)
```r
library(spatstat.geom); library(spatstat.core)

# 1. Load crash points (replace with your own shapefile)
pp  <- as.ppp(crash_coords, W = as.owin(study_window))

# 2. Estimate K with isotropic correction
Khat <- Kest(pp, correction = "iso")

# 3. CSR envelope with 99 simulations
env  <- envelope(pp, Kest, nsim = 99, correction = "iso")

# 4. Plot
plot(env, . - theo ~ r, legend = FALSE,
     main = "L(r) - r under CSR",
     ylab = "L(r) - r")
abline(h = 0, lty = 2)
```
### 6. Key Takeaways
CSR = homogeneous Poisson. Departures indicate spatial dependence or inhomogeneity.

Always edge‑correct; otherwise small windows falsely suggest inhibition.

Testing CSR is diagnostic, not a final model. If rejected, next fit an inhomogeneous Poisson or an interaction model.

# Lesson 3 – Intensity Estimation and the Inhomogeneous Poisson Model  
*A deep‑dive into first‑order structure*

---

## 1 Conceptual Foundation

### 1.1 Mean Measure vs. Intensity  
For a planar point process $\mathbf X$ observed in window $W\subset\mathbb R^{2}$,

$$
\Lambda(B) \;=\; \mathbb E[N(B)], \quad  
\lambda(s) \;=\; \frac{\partial^{2}\Lambda(B)}{\partial x\,\partial y}
\;\;\biggr|_{s\in B}.
$$

* $\Lambda(B)$ – expected **count** inside any Borel set $B$.  
* $\lambda(s)$ – expected **density** of events at location $s$ (points / unit area).

If events occur on a linear network $L$ of total length $|L|$ the denominator is **unit length**, and for space–time data it is **unit area × unit time**.

---

## 2 Non‑parametric Intensity: Kernel Smoothing

### 2.1 Fixed‑Bandwidth Estimator  

$$
\widehat{\lambda}_h(s) =
\frac{1}{h^2}
\sum_{i=1}^{N}
\frac{1}{e_i(s)}\,
K\left( \frac{ \lVert s - s_i \rVert }{ h } \right)
$$



where  

* $K(\cdot)$ is a radially symmetric kernel with $\int_{\mathbb R^{2}} K(\lVert u\rVert) \,du = 1$.  
* $h$ is the bandwidth (smoothing radius).  
* $e_i(s)=\int_{W} K_h(s-u)\,du$ is the **edge correction** weight.

> **Common kernels**  
> $K(r)=\frac{1}{2\pi}\exp(-\tfrac12 r^{2})$ (Gaussian) | $K(r)=\frac{3}{\pi}(1-r^{2})_{+}$ (Epanechnikov).

---

### 2.2 Bandwidth Selection

| Approach                     | Bandwidth‑selection criterion                                                                                      | Practical notes                                                                                                       | Key reference |
|------------------------------|--------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------|---------------|
| Rule‑of‑thumb                | $h_{\text{ROT}} = 0.15 \sqrt{|W|/N}$                                                                                | Extremely fast; assumes a homogeneous Poisson process, so it ignores spatial clustering.                              | (Diggle, 1985) |
| Likelihood cross‑validation  | $\displaystyle \max_{h}\sum_{i=1}^{N}\log\!\left[\widehat{\lambda}_{h,-i}(s_i)\right]$                              | Default in `spatstat::bw.ppl`; robust to inhomogeneity but computationally expensive for large $N$.                   | (Baddeley & Turner, 2005) |
| Least‑squares cross‑valid.   | $\displaystyle \min_{h} \int_{W} \left[ \widehat{\lambda}_{h}(s) - \widehat{\lambda}_{\text{pilot}}(s) \right]^2 \, ds$ | Sensitive to the pilot bandwidth; often undersmooths with sparse data.                                                | (Cronie & van Lieshout, 2016) |
| Plug‑in (variance matching)  | $\displaystyle h_{\text{PI}} = \arg\min_{h} \left| \widehat{\sigma}^{2}(h) - \sigma^{2}_{\text{target}} \right|$         | Yields stable estimates when $N$ is large; requires a pilot step to estimate the underlying Gaussian‑field variance.   | (Baddeley et al., 2016) |


> **Bias–variance trade‑off**  
> Small $h$ yields noisy “spikes”; large $h$ oversmooths and hides structure. Report both *map* and *bandwidth* in papers.

---

### 2.3 Adaptive Kernel Smoothing  

Set a pilot intensity $\tilde\lambda(s)$ with fixed $h_0$, then let  

$$h(s) = h_0 \biggl[\frac{\tilde\lambda(s)}{\text{median}\{\tilde\lambda(s)\}}\biggr]^{-1/2}.$$

Dense areas get a smaller bandwidth, preserving sharp hotspots; sparse areas get wider bandwidth, reducing variance.

---

## 3 Parametric Intensity: Log‑Linear Poisson Model

### 3.1 Model Specification  

$$
\lambda(s) \;=\; \exp\!\bigl(\beta_0 + \beta^{\top}Z(s)\bigr),
$$

where $Z(s) = \bigl(Z_1(s),\dots,Z_p(s)\bigr)$ are spatial covariate fields (elevation, population, AADT, land‑use fractions, …).  

* Interpretation: $\exp(\beta_k)$ is the multiplicative change in intensity for a one‑unit increase in $Z_k$ (holding others fixed).  
* Add an **offset** $\log(E(s))$ if observation effort $E(s)$ varies across $W$.

---

### 3.2 Likelihood via Berman–Turner Quadrature

Create $n_D$ dummy (integration) points $\{u_j\}$ with weights $w_j$.  
Let $y_i=1$ for events, $y_j=0$ for dummies.

$$\ell(\beta) \;=\;
\sum_{i=1}^{N}\beta^{\top}Z(s_i)
\;-\;
\sum_{j=1}^{n_D} w_j\,\exp\!\bigl[\beta^{\top}Z(u_j)\bigr].$$

This is equivalent to a weighted Poisson GLM:

```r
library(spatstat.geom); library(spatstat.core)
X  <- ppm(points_ppp,
          ~ elev + aadt + factor(landuse),
          covariates = list(elev = elev_raster,
                            aadt = aadt_raster,
                            landuse = landuse_raster))
summary(X)
```



