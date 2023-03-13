# Intrinsic spherical smoothing method based on generalized Bézier curves and sparsity inducing penalization
The method aims to smooth the data points observed along time points.
We use the shperical B\'ezier curves and local penalization scheme.
The penalization scheme induces the order selection and select an appropricate smoothness of spherical B \'e zier curve.

## Data

Let $\{ (t_n, y_n) \}_{n=1}^N$ be a set of data, 
where $t_n \in \mathbb{I}$ are the given time points 
and $y_n$ are the associated data points on $\mathbb{S}$. 
Our goal is to find a spherical B \'e zier curve that fits the data well.
