# Intrinsic spherical smoothing method based on generalized Bézier curves and sparsity inducing penalization
The method aims to smooth the data points observed along time points.
We use the shperical B\'ezier curves and local penalization scheme.
The penalization scheme induces the order selection and select an appropricate smoothness of spherical B \'e zier curve.

## Data

Let $\{ (t_n, y_n) \}_{n=1}^N$ be a set of data, 
where $t_n \in \mathbb{I}$ are the given time points 
and $y_n$ are the associated data points on $\mathbb{S}$. 
Our goal is to find a spherical B \'e zier curve that fits the data well.

## Spherical B \'e zier curve

A linear B\'ezier curve $\gamma_1: \mathbb{I} \to \mathbb{S}$ is simply a geodesic segment between $v$ and $w$ defined as

$$
\gamma_1(t; v, w) 
= \frac{\sin((1 - t) \theta)}{\sin\theta} v + \frac{\sin(t \theta)}{\sin\theta} w,
\quad t \in \mathbb{I},
$$

where $\theta = \mathsf{d}(v, w) \triangleq \arccos(v^\top w)$ denotes the great circle distance between $v$ and $w$.

The curve 

$$
\gamma_J (t; \xi) \triangleq \gamma_J(t; \xi_0, \ldots, \xi_J) \quad \text{for} \quad t \in \mathbb{I}
$$

is called the spherical B\'ezier curve of degree $J$ (order $J+1$) with control points $\xi_0, \ldots, \xi_J$.

## Objective function

The objective function to minimize is given by

$$
\ell^\lambda(\xi) = \ell(\xi) + \lambda \sum_{j=1}^{J-1} 
\left| D_j(\xi) \right| \for D_j(\xi) = \dot{\gamma}_1(0; \xi_j, \xi_{j+1}) - \dot{\gamma}_1(1; \xi_{j-1}, \xi_j) \quad \text{for} \quad \xi \in \Omega,
$$

where $\lambda > 0$ is the complexity parameter.

## Estimator 
For a fixed $\lambda$, the penalized intrinsic spherical B\'ezier (PISB) curve is defined as

$$
\gamma_J (\cdot; \hat{\xi}_\lambda),
$$

where

$$
\hat{\xi}_\lambda = \argmin_{\xi \in \Omega} \ell^\lambda(\xi). 
$$
