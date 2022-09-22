# The N-body Gravitational Problem

The study of planets and other objects in the solar system
was very important in the development of physics and Newton's
Law of Universal Gravitation.  
The gravitational two-body problem can be solved analytically.  
For bound systems, the solutions are elliptical orbits, 
as discovered by Kepler, while unbound systems have 
hyperboloidal orbits.
The $N$-body problem (for $N>2$) can not be solved analytically 
in general, and these systems can have very compicated--and
even chaotic--solutions.   This project explores some
examples of these solutions.

## Introduction

The force between two bodies with mass ``m_1`` and ``m_2`` has
a magnitude
```math
F = \frac{G m_1 m_2}{r^2},
```
where ``G`` is Newton's Gravitation constant, and ``r`` is the
distance between the two masses.  Using the position vectors for 
the two bodies, ``\mathbf{r}_1`` and ``\mathbf{r}_2``, we can
write this in vector form as
```math
\mathbf{F}_{12} = - \frac{G m_1 m_2}{|\mathbf{r}_1 - \mathbf{r}_2|^2}\hat{\mathbf{r}}_{12},
```
where ``\mathbf{F}_{12}`` is the force on body 1 due to the body 2.
and ``\hat{\mathbf{r}}_{12}`` is the unit vector from body 1 to 2
```math
\hat{\mathbf{r}}_{12} = \frac{\mathbf{r}_1 - \mathbf{r}_2}{|\mathbf{r}_1 - \mathbf{r}_2|}.
```
Newton's second law gives the equations of motion for two bodies as
```math
\begin{aligned}
\mathbf{F}_{\rm net} = m_1 \ddot{\mathbf{r}}_1 &=  \mathbf{F}_{12} 
= -\frac{G m_1 m_2}{|\mathbf{r}_1 - \mathbf{r}_2|^3}(\mathbf{r}_1 - \mathbf{r}_2)\\
m_2 \ddot{\mathbf{r}}_2 &=  \mathbf{F}_{21} 
= -\frac{G m_1 m_2}{|\mathbf{r}_2 - \mathbf{r}_1|^3}(\mathbf{r}_2 - \mathbf{r}_1)
\end{aligned}
```
By Newton's Third Law, ``\mathbf{F}_{21} = - \mathbf{F}_{12}``, as can be 
seen.

To generalize these equations for ``N`` bodies, we simply sum all of the
forces.  The force on body at position ``\mathbf{r}_a`` with mass ``m_a``
is
```math
\mathrm{F}_{\rm net} = m_a \ddot{\mathbf{r}}_a 
=  \sum_{b=1, b\neq a}^{\infty} 
-\frac{G m_a m_b}{|\mathbf{r}_a - \mathbf{r}_b|^3}(\mathbf{r}_a - \mathbf{r}_b)
```
As mentioned before, Hamiltonian systems have several advantages for
numerical work, such as numerical methods can conserve the energy and
angular momentum in problems like these.  The fundamental variables
of Hamiltonian systems are the coordinates, usually labeled with ``q``'s
and the canonical momentae, labeled with ``p``'s.  In Cartesian coordinates,
the canonical momentae are
```math
p_a = m_i\dot{x}_a.
```
Using this definition of the canonical momentum and Newton's second law
in the general form
```math
\frac{d\mathbf{p}}{dt} = \mathbf{F}_{\rm net},
```
we can write the equations of motion for this system as
```math
\begin{aligned}
\dot{\mathbf{q}}_a &= \frac{\mathbf{p}_a}{m_a}\\
\dot{\mathbf{p}}_a &= \sum_{b=1, b\neq a}^{\infty} 
-\frac{G m_a m_b}{|\mathbf{r}_a - \mathbf{r}_b|^3}(\mathbf{r}_a - \mathbf{r}_b).
\end{aligned}
```

## Numerical Solution

We can solve the ODEs using the solvers in the `DifferentialEquations.jl`
package.  We could use symplectic methods that conserve the total energy
and the angular momentum.  However, these methods require a constant time
step. Some three-body problems have very different dynamical time
scales at different times or for different bodies.  In this case, it may
be more efficient to use a high-order, adaptive integrator.

## Examples



## Extension to General Relativity

The Einstein equations of general relativity are complicated, nonlinear
partial differential equations.  However, in a weak gravitational field,
we can use approximations to simplify the equations.  One such method
for simplifying the equations is the __Post-Newton Method.__  
The post-Newtonian method is useful when the gravitational potential
is small *and* the particle speeds are also small compared to the
speed of light
```math
\left| \Phi\right|  = \left|\frac{Gm}{rc^2} \right| \ll 1, \qquad {\rm and}
\qquad \left(\frac{v^2}{c^2}\right) \ll 1.
```

Before introducing the post-Newtonian equations, we first begin with 
the Hamiltonian for Newtonian gravity.  Consider three bodies with 
mass ``m_a``, labelled with the index ``a = 1, 2, 3``.  The Hamiltonian
for this system is simply the total energy, expressed in terms of the 
positions of the particles ``\mathbf{r}_a`` and their momentae ``\mathbf{p}_a``.
We define the vector ``r_{ab} = \mathbf{r}_a - \mathbf{r}_b`` 
with magnitude ``r_{ab} = |\mathbf{r}_a - \mathbf{r}_b|,`` and the unit
vector ``\mathbf{n}_{ab} = \mathbf{r}_{ab}/r_{ab}``.
The square of the momentum is ``p_a^2 = \mathbf{p}_a \cdot \mathbf{p}_a``.
Hamilton's equations are
```math
\dot{\mathbf{q}}_a = \frac{\partial H}{\partial \mathbf{p}_a},
\qquad
\dot{\mathbf{p}}_a = -\frac{\partial H}{\partial \mathbf{q}_a}.
```

```math
H_N = \frac{1}{2}\sum_a \frac{p_a^2}{m_a} 
         - \frac{1}{2} \sum_{a,b\neq a} \frac{m_a m_b}{r_{ab}}
```

The post-Newtonian equations introduce corrections to the Newtonian
equations in an expansion in powers of ``\epsilon \sim \Phi \sim v^2/c^2``
```math
H = H_N + H_{1PN} + H_{2PN} + \cdots
```
The first-order corrections for ``N`` bodies is
```math
\begin{aligned}
H_{1PN} &= -\frac{1}{8} \sum_a m_a \left(\frac{p_a^2}{m_a^2}\right)^2
          - \frac{1}{4}\sum_{a,b\neq a} \frac{m_a m_b}{r_{ab}} 
          \left\{ 
             6\frac{p_a^2}{m_a^2} 
             - 7\frac{\mathbf{p}_a\cdot \mathbf{p}_b}{m_a m_b} 
          - \frac{(\mathbf{n}_{ab}\cdot\mathbf{p}_a)(\mathbf{n}_{ab}\cdot\mathbf{p}_b)}{m_a m_b}
          \right\}\\
       &\quad  + \frac{1}{2}\sum_{a, b\neq c, c\neq a}\frac{m_a m_b m_c}{r_{ab}r_{ac}}.
\end{aligned}
```
