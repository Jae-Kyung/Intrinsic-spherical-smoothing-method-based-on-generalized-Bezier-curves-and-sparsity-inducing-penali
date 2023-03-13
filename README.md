# Intrinsic spherical smoothing method based on generalized Bézier curves and sparsity inducing penalization
The method aims to smooth the data points observed along time points.
We use the shperical B\'ezier curves and local penalization scheme.
The penalization scheme induces the order selection and select an appropricate smoothness of spherical B \'e zier curve.

## Data

Let $\{ (t_n, y_n) \}_{n=1}^N$ be a set of data, 
where $t_n \in \mathbb{I}$ are the given time points 
and $y_n$ are the associated data points on $\mathbb{S}$. 
Our goal is to find a spherical B \'e zier curve that fits the data well.

## Spherical B\'ezier curve

A linear B\'ezier curve $\gamma_1: \mathbb{I} \to \mathbb{S}$ is simply a geodesic segment between $v$ and $w$ defined as

$$
\gamma_1(t; v, w) 
= \frac{\sin((1 - t) \theta)}{\sin\theta} v + \frac{\sin(t \theta)}{\sin\theta} w,
\quad t \in \mathbb{I},
$$

where $\theta = \mathsf{d}(v, w) \triangleq \arccos(v^\top w)$ denotes the great circle distance between $v$ and $w$.

We define

$$
\gamma_j(t; \xi_i, \ldots, \xi_{i+j}) =
\gamma_1 { t; \gamma_{j-1}(t; \xi_{i : i+j-1}), \gamma_{j-1}(t; \xi_{i+1 : i + j}) } \for t \in \mathbb{I},
$$

where $j = 2, \ldots, J$ and $i = 0,\ldots, J - j$. 

The curve 

$$
\gamma_J (t; \xi) \triangleq \gamma_J(t; \xi_0, \ldots, \xi_J) \quad \text{for} \quad t \in \mathbb{I}
$$

is called the spherical B\'ezier curve of degree $J$ (order $J+1$) with control points $\xi_0, \ldots, \xi_J$.

## Objective function

The objective function to minimize is given by

$$
\ell^\lambda(\xi) = 
\frac{1}{2} \sum_{n = 1}^N \mathsf{d}^2 (y_n, \gamma_J(t_n; \xi)) + \lambda \sum_{j=1}^{J-1} \left| D_j(\xi) \right| \quad \text{for} \quad D_j(\xi) = \dot{\gamma}_1(0; \xi_j, \xi_{j+1}) - \dot{\gamma}_1(1; \xi_{j-1}, \xi_j) \quad \text{for} \quad \xi \in \Omega,
$$

where $\lambda > 0$ is the complexity parameter.
As the complexity parameter $\lambda$ increases, the difference in the velocity vectors of the two geodesic segments becomes smaller.

## Estimator 
For a fixed $\lambda$, the penalized intrinsic spherical B\'ezier (PISB) curve is defined as

$$
\gamma_J (\cdot; \hat{\xi}_\lambda),
$$

where

$$
\hat{\xi}_\lambda = \text{argmin}_{\xi \in \Omega} \ell^\lambda(\xi). 
$$

## Figures
- Quadratic B\'ezier curve

![Fig1-2](https://user-images.githubusercontent.com/84615460/224720400-3556ce28-5948-461f-8875-a30ce0112cca.png)

- Penalization effects

The quadratic B\'ezier curve is to become the below B\'ezier curve when the complexity parameter increases.
Note that the final is the linear B\'ezier curve

![Fig2-2](https://user-images.githubusercontent.com/84615460/224721235-3e4a890b-2d66-4be6-a8ad-05a97794a37d.png)
![Fig2-3](https://user-images.githubusercontent.com/84615460/224721284-3d43f64f-b01b-4857-a544-d10abda48cf9.png)


