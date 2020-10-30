# RayTracingAngularQuadrature

## Description

Rays can be used to represent anisotropic transport of along straight lines (e.g. transport and attenuation of
radiation). Coupled applications often do not need the angular information propagated along each ray, but the
integral of the angular information over a solid angle $\Delta \hat{\Omega}$.  

The [modules/ray_tracing/index.md] module provides angular quadratures for computing integrals over solid angles.
These quadratures approximate the angular integral by a summation:

\begin{equation}
  \int\limits_{\Delta \hat{\Omega}} f (\hat{\Omega})d \hat{\Omega} \approx \sum\limits_{j=1}^J w_j f (\hat{\Omega}_j),
\end{equation}

where $f (\hat{\Omega})$ is the quantity that is propagated/computed along the ray, $\hat{\Omega}_j$ is angular direction of
this ray, and $w_j$ is an angular weight. The angular quadrature provides two essential functions: (1) it specifies which
directions $\hat{\Omega}_j$ to select for rays in order to propagate the quantity $f (\hat{\Omega})$, and (2) it provides
weights $w_j$ that allow approximating angular integrals of $f$.

The angular integral can be expanded into an integral over azimuthal angle $\phi$ and polar angle $\theta$. The polar angle is
measured between the direction $\hat{\Omega}$ and the z-axis, while the azimuthal angle is the angle between the projection
of the angular direction into the x-y plane and the positive x-axis. The infinitesimal element $d \hat{\Omega}$ is related to
$d\phi$ and $d\theta$ by:

\begin{equation}
d \hat{\Omega} = \sin \theta d\theta d\phi.
\end{equation}

The angular integral is often rewritten in terms of $\mu = \cos \theta$:

\begin{equation}
 \int_0^{2 \pi} d \phi \int\limits_{-\pi/2}^{\pi/2} sin \theta d \theta =
 \int_0^{2 \pi} d \phi \int\limits_{-1}^{1} d \mu.
\end{equation}

### Rotating the angular quadrature

It is convenient to create a quadrature rule that is oriented with respect to an arbitrary vector $\vec{n}$ instead of
the $z-axis$. We denote angle $\hat{\Omega}_j$ as angular direction with respect to the z-axis and $\hat{\Omega}_{j, \vec{n}}$
as angular direction with respect to the vector $\vec{n}$. We relate $\hat{\Omega}_{j, \vec{n}}$ and $\hat{\Omega}_j$ by:

\begin{equation}
   \hat{\Omega}_{j, \vec{n}} = \underline{R} \hat{\Omega}_j,
\end{equation}

where $\underline{R}$ is a rotation matrix given by:

\begin{equation}
\underline{R} = \left [
\begin{matrix}
n_x & -n_y & 0\\
n_y & n_x & 0 \\
0 & 0 & 1
\end{matrix}
\right ],
\end{equation}

for 2D geometries. For 3D geometries, we define the vectors $\vec{t}_x$ and $\vec{t}_y$, where $\vec{t}_x$
is chosen so that $\vec{t}_x \cdot \vec{n} = 0$ and $\vec{t}_y = \vec{n} \times \vec{t}_x$.  Then:

\begin{equation}
\underline{R} = \left [\vec{t}_x, \vec{t}_y, \vec{n} \right ].
\end{equation}

### Legendre-Chebyshev Quadrature

The Legendre-Chebyshev quadrature is a product quadrature. The $\theta$ direction is integrated using a Gauss-Legendre quadrature,
while the azimuthal direction is integrated using a Chebyshev quadrature. For an integral within a cone of width $\Delta \hat{\Omega}$, we have:

\begin{equation}\label{eq:cone_integral}
    \int_0^{2 \pi}  \int\limits_{1 - \Delta \mu}^{1} f(\phi, \mu) ~d \mu d \phi \approx \sum\limits_{i=1}^I  \sum\limits_{j=1}^J w_{C,i} w_{L} f(\phi_i, \mu_j),
\end{equation}

where

 1. Chebyshev weights: $w_{C,i} = 2 \pi / I$.

 2. Chebychev abscissae: $\phi_i = \frac{(2 i - 1) \pi}{I}$

 3. $\mu_j$ and $w_{L,j}$ are Gauss-Legendre abscissae and weights transformed to the integral over the range $[1 - \Delta \mu, 1]$.

### Angular Quadrature in 3D

In three-dimensional geometries, the quadrature is given by:

\begin{equation}
  \big \{\phi_i = \frac{(2 i - 1) \pi}{I}, w_{C,i} = 2 \pi / I, \mu_j, w_{L,j}, \hat{\Omega}_j = (\sqrt{1 - \mu_j^2} \cos \phi_i,\sqrt{1 - \mu_j^2} \sin \phi_i, \mu_j) \big \}_{i=1,..,I; j=1,..,J}
\end{equation}

The quadrature weights sum to:

\begin{equation}
 \sum\limits_{i=1}^I  \sum\limits_{j=1}^J w_{C,i} w_{L} = 2 \pi \Delta \mu.
\end{equation}

A half-range quadrature is obtained for $\Delta mu = 1$ and a full-range quadrature is obtained for $\Delta \mu = 2$.

### Angular Quadrature in 2D

In two-dimensional geometries, rays have only one angular degree of freedom: the polar angle $\theta$.
