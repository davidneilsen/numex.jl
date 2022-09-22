# Binary Orbits

Historically, the gravitational two-body problem was very important
in the development of physics. People have tracked the motion of planets 
across the night sky for thousands of years.  Predicting this 
planetary motion has been a concern of science since the time of 
ancient Greece, and it led eventually to the work of Copernicus, Kepler, 
Newton, and Einstein.

The mass of the our sun dominates
the mass of all the other planets in our solar system, so to a good
approximation planetary orbits can be found by simply considering
only the gravitational force of the sun.

## Newtonian Gravity

### Solution with Newton's Second Law

The gravitational force between two bodies with mass ``m_1`` and ``m_2`` has
a magnitude
```math
F = \frac{G m_1 m_2}{r^2},
```
where ``G`` is Newton's Gravitation constant, and ``r`` is the
distance between the two masses.  Let the sun be the first body,
``m_1 = m_{\odot},`` located at the origin, and the second body be
a planet with mass ``m_2 \ll m_{\rm sun},`` located a distance
``r`` from the origin.  

To solve the equations of motion given by Newton's second law, we need
to choose a coordinate system. In advanced mechanics courses, this
problem is solved analytically in polar coordinates. For numerical
solutions, Cartesian coordinates work well.  To write the equations
in vector form, let ``\mathbf{r} = x\mathbf{i} + y\mathbf{j}.``  
Newton's second law is
```math
m_2\ddot{\mathbf{r}} = -\frac{G m_1 m_2}{r^2} \hat{\mathbf{r}},
```
where ``r = \sqrt{x^2 + y^2}`` and 
``\hat{\mathbf{r}}`` is a unit vector in the radial direction
```math
\begin{aligned}
\hat{\mathbf{r}} &= \mathbf{i}\,\cos\theta + \mathbf{j}\, \sin\theta\\
  &= \frac{x}{r}\mathbf{i} + \frac{y}{r}\mathbf{j}.
\end{aligned}
```
To solve these equations numerically, we write them in as components
```math
\begin{aligned}
m_2 \ddot{x} &= -\frac{G m_1 m_2}{r^2}\frac{x}{r}\\
m_2 \ddot{y} &= -\frac{G m_1 m_2}{r^2}\frac{y}{r}.
\end{aligned}
```
in first-order form.  Perhaps the obvious way to reduce these equations
to first-order form is to introduce the velocities ``v_x = \dot x`` 
and ``v_y = \dot y,`` giving
```math
\begin{aligned}
\dot{x} &= v_x\\
\dot{y} &= v_y\\
\dot{v}_x &= -\frac{G m_1 x}{r^3}\\
\dot{v}_y &= -\frac{G m_1 y}{r^3}.
\end{aligned}
```

### Lagrangian

### Hamiltonian

## Numerical solutions

The equations of motion in first order form can be solved using
ODE solvers in the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) 
package.  

However, before solving the equations, we need to use
a unit system appropriate for numerical work.  In SI units,
the mass of the sun is ``m_\odot = 2\times 10^{30}`` kg,
Newton's constant is ``G=6.67\times 10^{-11}`` N m``{}^2``/kg``{}^2``,
and the distance between the Earth and Sun is ``r\simeq 1.5\times 10^{11}`` m.
These numbers have drastically different orders of magnitude, which 
can introduce round-off errors in numerical calculations.  
If we set the astronomical unit to one length unit, the solar mass to one
mass unit, and use the year as the time unit, then ``G=4\pi^2.``

## Relativistic Corrections

General relativity reduces to Newton's law of universal gravitation
in very weak gravitational fields.  The force law with the relativistic
correction for binary systems is
```math
F = -\frac{G m_1 m_2}{r^2} - \frac{3GM\ell^2}{m_2 c^2 r^4}
```
where ``\ell`` is the angular momentum per unit mass.  The equations
of motion can then be written
```math
\begin{aligned}
\dot{x} &= v_x\\
\dot{y} &= v_y\\
\dot{v}_x &= -\frac{G m_1 x}{r^3} - \frac{3GM\ell^2 x}{c^2r^5}\\
\dot{v}_y &= -\frac{G m_1 y}{r^3} - \frac{3GM\ell^2 y}{c^2r^5}.
\end{aligned}
```
The additional force term causes the elliptical orbits to precess.  The
precession of Mercury's orbit was first measured in 1859 by Urbain Le Verrier.
Most of the precession can be explained by perturbations from other
planets, such as Jupiter.  However, these perturbations did not explain 
the total observed precession.
While developing general relativity, Albert Einstein realized that relativistic
effects also cause orbital precession, and general relativity explains
the missing part of the Mercury's precession.  This is one of the three
classical tests of general relativity.



## References

1. Wang, F. Y.-H., "Relativistic orbits with computer algebra," *Am. J. Phys.* __72__ 1040 (2004); [DOI](https://doi.org/10.1119/1.1645284)
2. Boyer, T. H., "Unfamiliar trajectories for a relativistic particle in a Kepler or Coulomb potential," *Am. J. Phys.* __72__ 992 (2004); [DOI](https://doi.org/10.1119/1.1737396)
3. Lemmon, T. J. and Mondragon, A. R., "Kepler's Orbits and Special Relativity in Introductory Classical Mechanics," arXiv:1012.5438[astro-ph.EP]; [DOI](https://doi.org/10.48550/arXiv.1012.5438)
