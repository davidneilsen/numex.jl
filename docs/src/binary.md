# Binary Orbits

People have tracked the motion of planets 
across the night sky for thousands of years.  Predicting this 
planetary motion has been a concern of science since the time of 
ancient Greece, and it led eventually to the work of Copernicus, Kepler, 
Newton, and Einstein.

In this project we begin with a simple model using Newtonian gravity.
We then consider corrections to this model from general relativity.
Similar models are used by LIGO to search for gravitational waves.

## Newtonian Gravity

Newton's universal law of gravitation gives the force between
two bodies with mass ``m_1`` and ``m_2`` has a magnitude
```math
F = \frac{G m_1 m_2}{r^2},
```
where ``G`` is Newton's Gravitation constant, and ``r`` is the
distance between the two masses.

### Approximate solution with Newton's Second Law

The mass of the our sun dominates the mass of all the other planets
in our solar system.  To a good approximation planetary orbits
can be found by simply considering the motion of a single planet
around a fixed point at the center of the sun.  This approximation
works well when the planet's mass, ``m``, is much less than the
mass of the sun, ``M_\odot,`` so that the motion of the sun can 
be ignored.


We first choose a coordinate system.  The sun is fixed at the origin,
and the planet's position is ``\mathbf{r}(t)``.
The orbit is confined to a plane, so we could use Cartesian coordinates, 
``(x,y)``, but polar coordinates
```math
r = \sqrt{x^2+ y^2}, \qquad \tan\theta = \frac{y}{x}.
```
are a natural choice as the gravitational potential is a function of ``r``.
The calculation of vector derivatives such as the velocity, however, are
easiest to perform in Cartesian coordinates
```math
x = r\cos\theta, \qquad y = r\sin\theta.
```
Finally, as a short-hand notation, we use a dot to represent a time derivative
```math
\dot{x} \equiv \frac{dx}{dt}.
```

To solve for the planet's motion, we could use Newton's second law directly.
However, conserved quantities often simplify the analysis, so let's see if
we can use them.  The total energy in vector form is
```math
E = \frac{1}{2} m \dot{\mathbf{r}}^2 - \frac{GMm}{r}.
```
The kinetic energy is a function of the velocity squared
```math
\dot{\mathbf{r}}^2 = \dot{x}^2 + \dot{y}^2.
```
We write this equation in terms of polar coordinates, and using some trig
identities, we get
```math
\dot{\mathbf{r}}^2 = \dot{r}^2 + r^2\dot{\theta}^2.
```
Then the energy can written in polar coordinates as
```math
E = \frac{1}{2} m \left(\dot{r}^2 + r^2\dot{\theta}^2\right) - \frac{GMm}{r}.
```

Energy is not the only conserved quantity in gravitational interactions.
The angular momentum is also conserved. The angular momentum is
```math
\mathbf{J} = \mathbf{r} \times \mathbf{p}.
```
The planet is confined to the plane with angular velocity ``r\dot{\theta},``
so the angular momentum only has a component normal to the plane.
```math
J = r \cdot(mr\dot{\theta}) = mr^2\dot{\theta}.
``` 

We can use the angular momentum to eliminate ``\dot{\theta}`` from the energy equation as ``\dot{\theta} = \frac{J_z}{mr^2},`` giving 
```math
E = \frac{1}{2} m \left(\dot{r}^2 + r^2 \dot{\theta}^2 \right) - \frac{GMm}{r},
```
The second and third terms on the RHS depend now only on position, and not
on velocity.  We combine these two terms into the effective potential
```math
V_{\rm eff} = \frac{J_z^2}{2m r^2} - \frac{GMm}{r},
```
giving 
```math
E = \frac{1}{2} m \dot{r}^2 + V_{\rm eff}(r).
```
Notice that ``\theta`` terms have now been eliminated from this equation
for the energy, making this now effectively a one-dimensional problem.

To find the equations of motion, we can differentiate the energy equation
with respect to time
```math
\frac{dE}{dt} = m \dot{r} \ddot{r} + \frac{dV_{\rm eff}(r)}{dr}\dot{r} = 0.
```
In writing this equation, we have used the chain rule.  For example,
```math
\frac{dV}{dt} = \frac{dV}{dr}\frac{dr}{dt} = \frac{dV}{dr}\dot{r}.
```
The derivative of the energy equation gives the equations of motion,
equivalent to Newton's second law
```math
m\ddot{r} = -\frac{dV_{\rm eff}(r)}{dr}.
```
We write these equations in first-order form, using the radial momentum
``p_r = m\dot{r}`` and the angular momentum ``J = mr^2\dot{\theta}``
```math
\begin{aligned}
\dot{r} &= \frac{p_r}{m}\\
\dot{\theta} &= \frac{J_z}{mr^2}\\
\dot{p}_r &= \frac{J_z^2}{m r^3} - \frac{GMm}{r^2}\\
\dot{J} &= 0.
\end{aligned}
```
These equations hold in the approximation that the first body does not
move, and that ``m_1\gg m_2``.

### Lagrangian

A more careful analysis of the binary problem for arbirary masses is easily
done with the Lagrangian formulation of mechanics.  The derivation
can be found in books on classical dynamics.  The Lagrangian is
found to be
```math
L = \frac{1}{2}\mu\dot{r}^2 + \frac{J_z}{2\mu r^2} + \frac{Gm_1 m_2}{r}.
```
where ``\mu`` is the reduced mass
```math
\mu = \frac{m_1 m_2}{m_1 + m_2}
```
and ``J`` is the conserved ``z``-component of the angular momentum
```math
J = \mu r^2 \dot{\phi}.
```


## Relativistic Corrections

### The Post-Newtonian Approximation

The Hamiltonian for the binary system is often derived using normalized
variables.  Working in the center-of-mass frame, 
``\mathbf{p}_1 + \mathbf{p}_2 = 0``, we define
```math
\mathbf{r} = \frac{\mathbf{x}_1 - \mathbf{x}_2}{GM},\qquad
\mathbf{p} = \frac{\mathbf{p_1}}{\mu} = -\frac{\mathbf{p_2}}{\mu}, \qquad
\hat{t} \equiv \frac{t}{GM}, \qquad
\hat{H} \equiv \frac{H}{\mu},
```
where
```math
M \equiv m_1 + m_2, \qquad \mu \equiv \frac{m_1 m_2}{M}, \qquad \nu \equiv \frac{\mu}{M}.
```
```math
\dot{\mathbf{r}} = \frac{d\mathbf{r}}{d\hat{t}} = \frac{d(\mathbf{x}_1 - \mathbf{x}_2)}{dt}, 
\qquad \dot{\mathbf{p}} = \frac{d\mathbf{p}}{d\hat{t}} = \frac{G}{\nu}\frac{d\mathbf{p}_1}{dt}.
```
```math
\hat{H}(\mathbf{r},\mathbf{p}) = \hat{H}_{N}(\mathbf{r},\mathbf{p})
    + \frac{1}{c^2}\hat{H}_{1PN}(\mathbf{r},\mathbf{p})
    + \frac{1}{c^4}\hat{H}_{2PN}(\mathbf{r},\mathbf{p}).
```
```math
\begin{aligned}
\hat{H}_{N}(\mathbf{r},\mathbf{p}) &= \frac{\mathbf{p}^2}{2} - \frac{1}{r},\\
\hat{H}_{1PN}(\mathbf{r},\mathbf{p}) &= \frac{1}{8}(3\nu -1)\left(\mathbf{p}^2\right)^2
   - \frac{1}{2}\left[ \left(3 + \nu \right)\mathbf{p}^2 
   + \nu(\mathbf{n}\cdot\mathbf{p})^2\right]\frac{1}{r} + \frac{1}{2r^2},\\
\hat{H}_{2PN}(\mathbf{r},\mathbf{p}) &= \frac{1}{16}\left(1 - 5\nu + 5\nu^2\right)\left(\mathbf{p}^2\right)^3 \\
 &\quad   + \frac{1}{8}\left[ \left(5 - 20\nu - 3\nu^2\right) \left(\mathbf{p}^2\right)^2  
         -2\nu^2(\mathbf{n}\cdot\mathbf{p})^2\mathbf{p}^2 
         - 3\nu^2(\mathbf{n}\cdot\mathbf{p})^4\right]\frac{1}{r}\\
&\quad + \frac{1}{2}\left[ \left(5 + 8\nu\right)\mathbf{p}^2  + 3\nu(\mathbf{n}\cdot\mathbf{p})^2 \right]
       \frac{1}{r^2} - \frac{1}{4}(1+3\nu)\frac{1}{r^3}.
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

## Numerical solutions

The equations of motion in first order form can be solved using
ODE solvers in the [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) 
package.  

The repository has an example Jupyter notebook [`PostNewtonianBinaryOrbits.ipynb`](https://github.com/davidneilsen/numex.jl/blob/main/examples/PostNewtonianBinaryOrbits.ipynb) that solves the binary
equations in Newtonian and post-Newtonian gravity.

## References

T. Damour, P. Jaranowski, and G. Schäfer, "Dynamical invariants for general relativistic two-body systems at the third post-Newtonian approximation," *Physical Review* **D62** 044024 (2000). [Journal](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.62.044024), or [arXiv](https://arxiv.org/abs/gr-qc/9912092).
