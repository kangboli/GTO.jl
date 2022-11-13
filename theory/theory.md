---
title: Gaussian Integrals
author: Kangbo
---


Much of this document follows [Boys 1949](https://royalsocietypublishing.org/doi/10.1098/rspa.1950.0036) and
[McMurchie & Davidson 1978](https://doi.org/10.1016/0021-9991\(78\)90092-X).

## The electrostatic potential of a charged sphere.

The most basic building block of the Gaussian integral is the electrostatic
potential of a charged sphere. In general, for a charge distribution
$\sigma \rho(\mathbf{r})$, where $\sigma$ is the total charge and $\rho(\mathbf{r})$
is a normalized density, the potential energy can be evaluated as
$$
V(\mathbf{r}) = \sigma \int \mathrm{d} \mathbf{r'} \frac{\rho(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|}
$$

For a uniformly charged sphere in 3D, the potential is $\rho(\mathbf{r}) =
\sigma \delta(r-R)$ its potential at $\mathbf{r}$ can be analytically evaluated
in the spherical coordinate.

![Spherical Charge](figures/gto.svg)

$$
\begin{aligned}
V(r)
&= \sigma \int \mathrm{d}r' r'^2\sin \theta \mathrm{d} \theta \mathrm{d} \phi  \frac{\delta(r' - R)} {R^2 + r^2 - 2 Rr \cos \theta}
= \sigma \int \frac{R^2 \sin \theta \mathrm{d} \theta \mathrm{d} \phi}{ R^2 + r^2 - 2 Rr \cos \theta} \\
&= \frac{2\pi \sigma}{r} \left(\sqrt{(R + r)^2} - \sqrt{(R -r)^2} \right)
=\begin{cases}
\frac{R^2 \sigma}{r} \quad \text{if }R < r\\
R \sigma \quad \text{if } R > r
\end{cases}
\end{aligned}
$$

## Gaussian as stacked concentric spheres.

A isotropic Gaussian distribution in 3D can now be viewed as many concentric spheres stacked together
$$
\rho(\mathbf{r}) = \exp(-\alpha r^2) = \int \mathrm{d} R \exp(-\alpha R^2) \delta(r-R).
$$
The potential is then
$$
\begin{aligned}
V(\mathbf{r}) &= 4\pi \int_0^r \mathrm{d} R \exp(-\alpha R^2) \frac{\sigma R^2}{r}
+ 4\pi \int_r^{\infty} \mathrm{d} R \exp(-\alpha R^2) \sigma R \\
&= \frac{2 \pi}{\alpha r} \int_{0}^r \exp(-\alpha R^2) \mathrm{d} R
= \frac{2 \pi}{\alpha^{3/2} r} \int_0^{\sqrt{\alpha} r} e^{-t^2} \mathrm{d}t
= \frac{\pi^{3/2} \mathrm{erf}\left(\sqrt{\alpha } r\right)}{2 \alpha ^{3/2} r}.
\end{aligned}
$$
Writing the integral as an error functions is intuitive, but 
a different form is more frequently used with $u = t/(\sqrt{\alpha} r)$.
The rationale for this form appears to be making the derivatives easier.
$$
V(r) = \frac{2 \pi}{\alpha} \int_0^1 e^{-\alpha r^2 u^2}\mathrm{d} u \triangleq
\frac{2 \pi}{\alpha} F_0(\alpha r^2). 
$$

The nuclear integral can be treated as a point charge in a Gaussian potential, and 
the above integral should suffice (without the angular part). For electron-electron 
integrals, we would need the potential energy of Gaussian in a Gaussian potential.

## A Gaussian in a Gaussian potential 

Given two Gaussians $\rho_1(\mathbf{r})$ and $\rho_2(\mathbf{r})$ centered at $\mathbf{r}_1$
and $\mathbf{r}_2$, the Coulomb energy 
in between is 
$$
\int \mathrm{d} \mathbf{r} \mathrm{d} \mathbf{r}'
\rho_1(\mathbf{r}) \frac{1}{|\mathbf{r} - \mathbf{r}'|} \rho_2(\mathbf{r}') 
$$
This integral is performed in the two center bipolar coordinate. This
coordinate is not very common where the coordinate of a point is specified by
its distances to $\mathbf{r}_1$ and $\mathbf{r}_2$, which we denote $p$ and
$q$, and an azimuthal angle $\phi$. The volume element is the shaded region,
whose area is $\mathrm{d}p \mathrm{d}q/\sin(\theta)$ with an thickness of 
$h \mathrm{d}\phi$. $h$ can be expressed through the area of the triangle as 
$h = pq \sin(\theta) / a$. Therefore, the volume element turns out to be 
$(p q/a) \mathrm{d}p \mathrm{d}q \mathrm{d}\phi$.


![Volume Element](figures/volume_elem.svg)

Perform the integral over in the bipolar coordinates gives
$$
\begin{aligned}
J &= \int \mathrm{d} \mathbf{r} \rho_1(\mathbf{r}) V(\mathbf{r})\\
&= \int (p q/a) \mathrm{d}p \mathrm{d}q \mathrm{d} \phi \rho_1(\mathbf{r})
\frac{2 \pi}{\alpha_1 |\mathbf{r} - \mathbf{r}_2|} \int_{0}^{|\mathbf{r} - \mathbf{r}_2|} \exp(-\alpha_1 R^2) \mathrm{d} R\\
&= \frac{(2 \pi)^2}{\alpha_1 a}  \int p \exp(-\alpha_2 p^2) \mathrm{d}p \mathrm{d}q  
\int_{0}^{q} \exp(-\alpha R^2) \mathrm{d} R\\
&= \frac{(2 \pi)^2}{\alpha_1 a}  
\left(\int_0^a \mathrm{d}q \int_{a-q}^{a+q} \mathrm{d}p + 
\int_a^{\infty} \mathrm{d}q \int_{q-a}^{q+a} \mathrm{d}p\right)
p \exp(-\alpha_2 p^2) 
\int_{0}^{q} \exp(-\alpha R^2) \mathrm{d} R\\
&= \frac{(2 \pi)^2}{2 \alpha_1 \alpha_2 a}
\left(\int_0^a \mathrm{d}q + \int_a^{\infty} \mathrm{d}q\right)
(\exp(-\alpha_2  (a-q)^2)-\exp(-\alpha_2 (a+q)^2))
\int_{0}^{q} \exp(-\alpha R^2) \mathrm{d} R\\
&= \frac{(2 \pi)^2}{2 \alpha_1 \alpha_2 a}
\int_{-\infty}^{\infty} \mathrm{d}q \exp(-\alpha_2  (a-q)^2)
\int_{0}^{q} \exp(-\alpha R^2) \mathrm{d} R\\
\end{aligned}
$$
To continue, we use a trick that has never been intuitive to me
$$
\begin{aligned}
J \frac{2 \alpha_1 \alpha_2 a}{(2\pi)2} =&\int_{-\infty}^{a} \mathrm{d} a' \int_{-\infty}^{\infty} \mathrm{d}q \frac{\partial }{\partial a'}\exp(-\alpha_2  (a'-q)^2)
\int_{0}^{q} \exp(-\alpha_1 R^2) \mathrm{d} R\\
=&\int_{-\infty}^{a} \mathrm{d} a' \int_{-\infty}^{\infty} \mathrm{d}q \frac{\partial }{\partial q}\exp(-\alpha_2  (a'-q)^2)
\int_{0}^{q} \exp(-\alpha_1 R^2) \mathrm{d} R\\
=&\int_{-\infty}^{a} \mathrm{d} a' 
\left(
\left.\exp(-\alpha_2  (a'-q)^2) \int_{0}^{q} \exp(-\alpha_1 R^2) \mathrm{d} R \right|_{-\infty}^{\infty} +\right.\\
& \left.\int_{-\infty}^{\infty} \mathrm{d}q \exp(-\alpha_2  (a'-q)^2)
\exp(-\alpha_1 q^2) \right)\\
=&\int_{-\infty}^{a} \mathrm{d} a' 
\sqrt{\frac{\pi}{\alpha_1 + \alpha_2}} \exp\left(-\frac{\alpha_1 \alpha_2 a'^2}{\alpha_1 + \alpha_2}\right)
\end{aligned}
$$

Write this in terms of $F_0$ yields
$$
J = \frac{2 \pi^{5/2}}{\alpha_1 \alpha_2 \sqrt{\alpha_1 + \alpha_2}} F_0\left(
\frac{a^2 \alpha_1 \alpha_2}{\alpha_1 + \alpha_2}
\right)
$$


## Gaussian with a polynomial

The derivative of a Gaussian w.r.t. its center is a Gaussian with a polynomial.
$$
\frac{\partial }{\partial x_C}
\exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2) =2 \alpha (x - x_C) \exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2)
$$
More generally, the $n^{\mathrm{th}}$ gives the $n^{\mathrm{th}}$ Hermite polynomial $\alpha^{n/2} H_n(\sqrt{\alpha} (x -x_C))$ by the definition of Hermite polynomials
$$
\frac{\partial^n }{\partial x_C^n} \exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2) = \alpha^{j/2} H_n(\sqrt{\alpha} (x - x_C)) \exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2).
$$
The Hermite polynomials form a complete basis for the vector space of polynomials up to some order,
so we can expand any polynomial in terms of Hermite polynomials as
$$
\begin{aligned}
P(x, y, z) \exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2)
&= \sum_{i,j,k} c_{i,j,k} \alpha^{(i+j+k)/2} H_i(\sqrt{\alpha}(x-x_C))H_j(\sqrt{\alpha}(y-y_C))H_k(\sqrt{\alpha}(z-z_C)) \exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2) \\
&= \sum_{i,j,k} c_{i,j,k} \frac{\partial^i }{\partial x_c^i}
\frac{\partial^j }{\partial y_c^j}
\frac{\partial^k }{\partial z_c^k}
\exp(-\alpha (\mathbf{r} - \mathbf{r}_C)^2).
\end{aligned}
$$
The coefficients $c_{i,j,k}$ can be found by a matching the monomial
coefficients. This matching is just a basis conversion and can be solved by a
linear solve, which turns out to be a back-substitution for most polynomials
basis including the Hermite polynomial basis. The back-substitution is often
done implicitly through the recursive relations of the Hermite polynomials, which
is conceptually more obscure but may make a difference on performance.

The potential of a Gaussian (we assume Gaussians to be with a polynomials from here on)
is then
$$
V(\mathbf{r}) = \sum_{i,j,k} c_{i,j,k} \frac{\partial^i }{\partial x_c^i}
\frac{\partial^j }{\partial y_c^j}
\frac{\partial^k }{\partial z_c^k}
\frac{2\pi}{\alpha} F(\alpha (\mathbf{r} - \mathbf{r}_C)^2)
\triangleq \sum_{i,j,k} c_{i,j,k} R_{i,j,k}(\alpha (\mathbf{r} - \mathbf{r}_C)^2).
$$
Note that $R_{i,j,k}$ is a universal set of scalar functions that can potentially be 
fitted instead of analytically evaluated. That being said, $R_{i,j,k}$ is generally 
recursively evaluated for efficient block evaluation.
$$
\begin{aligned}
T &= \alpha (a^2 + b^2 + c^2),\\
R_{0,0,0,n}(T) &= (-2 \alpha)^n F_n(T),\\
F_n(T) &= \int_0^1 u^{2n} e^{-T u^2} \mathrm{d} u,\\
R_{0,0,k+1,n}(T) &= c R_{0,0,k, n+1}  + k R_{0,0,k-1,n+1},\\
R_{0,j+1,k,n}(T) &= b R_{0,j,k, n+1}  + j R_{0,j-1,k,n+1},\\
R_{i+1,j,k,n}(T) &= b R_{i,j,k, n+1}  + i R_{i-1,j,k,n+1},
\end{aligned}
$$
where $F_n(T)$ is to be numerically evaluated as the incomplete Gamma function.

The Coulomb energy between two Gaussians can be similarly derived as
$$
J = \frac{2 \pi^{5/2}}{\alpha_1 \alpha_2 \sqrt{\alpha_1 + \alpha_2}}
\sum_{i,j,k,i',j',k'} c_{i,j,k} d_{i',j',k'}
(-1)^{i'+j'+k'} R_{i + i', j+j', k+k'}\left(\frac{\left(\mathbf{r}_1 - \mathbf{r}_2\right)^2 \alpha_1 \alpha_2}{\alpha_1 + \alpha_2} \right)
$$


## Product of two (or more) Gaussians

The product of two Gaussians is a Gaussian, and the new coefficients can be found by completing the squares.
$$
\begin{aligned}
&P_1(x, y, z) \exp(-\alpha_1 (\mathbf{r - \mathbf{r}_1}))
P_2(x, y, z) \exp(-\alpha_2 (\mathbf{r - \mathbf{r}_2})) \\
=& P_1(x, y, z) P_2(x, y, z) 
\exp\left(-\frac{\alpha_1 \alpha_2}{\alpha_1 + \alpha_2} (\mathbf{r}_1 - \mathbf{r}_2)^2 \right)
\exp\left(-(\alpha_1 + \alpha_2) \left(\mathbf{r} - \frac{\alpha_1 \mathbf{r_1} + \alpha_2 \mathbf{r}}{\alpha_1 + \alpha_2} \right)\right)\\
\triangleq& P_3(x,y,z) \exp(-\alpha_3 (\mathbf{r} - \mathbf{r}_3)^2)\\
\end{aligned}
$$
