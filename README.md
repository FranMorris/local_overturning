# Local overturning code
This repository contains code which calculates the zonal and meridional components of overturning mass flux, defined by [Schwendike et al. (2014)](https://doi.org/10.1002/2013JD020742).

## Local overturning definitions
In pressure coordinates, the continuity equation is $\nabla^2_p\chi+\frac{\partial\omega}{\partial p}=0$, where $p$ is pressure, $\nabla_p$ is the gradient operator on a pressure surface, and $\chi$ is the velocity potential defining the divergent part of the horizontal wind: $\mathbf{u_\text{div}}=\nabla_p\chi$.

Let $\mu$ be a potential function such that $\chi=\frac{\partial \mu}{\partial p}$. Combining this with the continuity equation, we obtain the expression $\nabla^2_p\mu=-\omega$, a Poisson equation.

Let $\mathbf{\psi}$ be a vector streamfunction $\mathbf{\psi}=-\nabla_p\mu$. Then, $\nabla_p\cdot\mathbf{\psi}=\omega$, so $\mathbf{\psi}$  points towards regions of ascending motion. The divergent wind can be separated into its zonal and meridional componnnts as $\mathbf{u_\text{div}}=(u_\lambda,u_\phi)$ and from the previous equations can then be defined as 

$$\mathbf{u_\text{div}}=\nabla_p\chi=\nabla_p\frac{\partial \mu}{\partial p}=\frac{\partial}{\partial p}\left(\nabla_p\mu\right)=\frac{\partial}{\partial p}\left(-\mathbf{\psi}\right)=-\left(\frac{\partial \psi_\lambda}{\partial p},\frac{\partial \psi_\phi}{\partial p}\right).$$

So what we have is two components of this vector streamfunction: $\psi_\lambda$ and $\psi_\phi$, the zonal and meridional componetns of the overturning streamfunction, respectively.

In the zonal direction, the continuity equation is satisfied by $\frac{1}{R_E}\frac{\partial u_\lambda}{\partial \lambda}+\frac{\partial}{\partial p}\left(\omega_\lambda\cos\phi\right)=0$ and in the meridional direction this is $\frac{1}{R_E}\frac{\partial}{\partial \phi}\left(u_\phi\cos{\phi}\right)+\frac{\partial}{\partial p}\left(\omega_\phi\cos{\phi}\right)=0$. The vertical motion $\omega$ must be the sum of the partitioned vertical motion, $\omega=\omega_\lambda+\omega_\phi$ and from these we can define overturning mass flux in the zonal and meridional directions: $m_\lambda=-\frac{\omega_\lambda\cos{\phi}}{g}$ and $m_\phi=-\frac{\omega_\phi\cos{\phi}}{g}$ . 


## What this code does
The code in this repository does the following:
1. solves the Poisson equation $\nabla_p^2\mu=-\omega$ using a global $\omega$ field on a pressure surface to calculate the potential function $\mu$.
2. calculates zonal and meridional components of the overturning streamfunction, $\psi_\lambda$ and $\psi_\phi$ by takinggradients of $\mu$ on the sphere
3. calculates components of overturning vertical velocity are computed by taking further gradients of the $\psi$ components: $\omega_\lambda=\frac{\partial\psi_\lambda}{\partial\lambda}$ and $\omega_\phi=\frac{\partial\psi_\phi}{\partial\phi}$.
4. includes a function to calculate corresponding overturning mass flux components: $m_i=-\frac{\omega_i\cos{\phi}}{g}$.

It also features functionality to splice a local area model (LAM) domain into a global domain with a matching grid, so that the calculation can be carried out on regional domains. (Regridding to the same domain is left as an exercise for the user.) This should be done using the LAM's driving global model. There may be some harmonic artefacts that will emerge beyond the LAM domain when taking the difference between the solution for the global domain and the LAM spliced into the global domain.

### Pressure surfaces
The equations are set up on pressure surfaces, so the calculations are done for a single pressure level at a time. You can loop over pressure levels to calculate this for the full atmospheric column. (Exercise caution at near-surface levels which may be obscured by orography and interpolated underground.)

### Solving a Poisson equation on the sphere
$\nabla^2\mu=\omega$ is a Poisson equation.

Spherical harmonics $Y_{lm}$ are the eigenfunctions of the Laplacian on the sphere, with eigenvalues $\frac{-l(l+1)}{R_E^2}$. So in spherical coordinates, $\nabla_p^2 Y_{lm}= -\frac{l(l+1)}{R_E^2} Y_{lm}$. 

For our Poisson equation, we can say that if our $\mu$ can be decomposed into some linear combination with coefficients of the cosine and sine terms within the spherical harmonics, e.g. $\mu=\sum_{l,m}{\mu_{lm} Y_{lm}(\phi,\theta)}$, then these coefficients are simply a factor out of being the same as $\omega=-\frac{l(l+1)}{R_E^2} \sum_{l,m} \omega_{lm} Y_{lm}$. So if you can get the $\omega_{lm}$ coefficients, you can get the $\mu_{lm}$ coefficients quite easily.

`SHExpandDH` is part of pyshtools which takes a Fourier transform of a grid of data points, then a weighted Legendre sum, to produce the coefficients. We apply this to the `omega` field and then loop through values of $l$ and divide it by $-l(l+1)$ to get the coefficients for $\mu$ which can then be expanded into the full $\mu$ field with `pyshtools.SHCoeffs.from_array()`. From this we can take gradients to obtain the streamfunction and vertical velocity components.


### Truncation
By choosing to only expand to a certain value of $l$, we can artificially "truncate" the solution - only taking the harmonics up to a certain number. In this code, `compute_overturning()` allows you to specify a truncation, a list of truncations, or not specify a truncation (it will default to calculating to the highest possible truncation for the grid).

### Units
This formulation solves on a unit sphere. Therefore the streamfunctions $\psi$ that we have are in units of $\text{Pa}\,\text{s}^{-1}$ rather than $\text{Pa}\,\text{m}\,\text{s}^{-1}$: multiply them by $R_E$ to give them physical units. The $\omega$ components remain in $\text{Pa}\,\text{s}^{-1}$.
