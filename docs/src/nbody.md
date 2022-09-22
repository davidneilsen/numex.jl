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
the two bodies, ``\mathrm{r}_1`` and ``\mathrm{r}_2``, we can
write this in vector form as
```math
\mathbf{F}_{12} =  -{G m_1 m_2}{|\mathrm{r}_1 - \mathrm{r}_2|^2} 
\hat \mathrm{r}_{12},
```
where ``\mathrm{F}_{12}`` is the force on body 1 due to the body 2.
and ``\mathrm{r}_{12}`` is the unit vector from body 1 to 2
```math
\mathrm{r}_{12} = \frac{\mathrm{r}_1 - \mathrm{r}_2}{|\mathrm{r}_1 - \mathrm{r}_2|}.
```
The equations of mtion for two bodies are then
```math
\begin{aligned}
m_1 \ddot \mathrm{r}_1 &=  \mathrm{F}_{12} 
= -{G m_1 m_2}{|\mathrm{r}_1 - \mathrm{r}_2|^3}(\mathrm{r}_1 - \mathrm{r}_2)\\
m_2 \ddot \mathrm{r}_2 &=  \mathrm{F}_{21} 
= -{G m_1 m_2}{|\mathrm{r}_2 - \mathrm{r}_1|^3}(\mathrm{r}_2 - \mathrm{r}_1)
\end{aligned}
```
By Newton's Third Law, ``\mathrm{F}_{21} = - \mathrm{F}_{12}``, as can be 
seen.

To generalize these equations for ``N`` bodies, we simply sum all of the
forces.  The force on body at position ``\mathrm{r}_i`` with mass ``m_i``
is
```math
m_i \ddot \mathrm{r}_i 
=  \sum_{\substack{j=1 \\ j\neq i}}^{\infty} 
-{G m_i m_j}{|\mathrm{r}_i - \mathrm{r}_j|^3}(\mathrm{r}_i - \mathrm{r}_j)
```
As mentioned before, Hamiltonian systems have several advantages for
numerical work, such as numerical methods can conserve the energy and
angular momentum in problems like these.  The fundamental variables
of Hamiltonian systems are the coordinates, usually labeled with ``q``'s
and the canonical momentae, labeled with ``p``'s.  In Cartesian coordinates,
the canonical momentae are
```math
p_i = m_i\dot x_i.
```
Using this definition of the canonical momentum and Newton's second law
in the general form
```math
\frac{d\mathrm{p}}{dt} = \mathrm{F}_{\rm net},
```
we can write the equations of motion for this system as
```math
\begin{aligned}
\dot \mathrm{q}_i &= \frac{\mathrm{p}_i}{m_i}\\
\dot \mathrm{p}_i &= \sum_\limits{j=1 \\ j\neq i}^{\infty} 
-{G m_i m_j}{|\mathrm{r}_i - \mathrm{r}_j|^3}(\mathrm{r}_i - \mathrm{r}_j).
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
\left| \Phi  = \frac{Gm}{rc^2} \right| \ll 1 \mbox{and} \frac{v^2}{c^2} \ll 1.
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
\dot \mathbf{q}_a = \frac{\partial H}{\partial \mathbf{p}_a},
\qquad
\dot \mathbf{p}_a = -\frac{\partial H}{\partial \mathbf{q}_a}.
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
The first-order corrections for three bodies is
```math
H_{1PN} = -\frac{1}{8} \sum_a m_a \left(\frac{p_a^2}{m_a^2}\right)^2
          - \frac{1}{4}\sum_{a,b\neq a} \frac{m_a m_b}{r_{ab}} 
          \left\{ 
             6\frac{p_a^2}{m_a^2} 
             - 7\frac{\mathbf{p}_a\cdot \mathbf{p}_b}{m_a m_b} 
          - \frac{(\mathbf{n}_{ab}\cdot\mathbf{p}_a)(\mathbf{n}_{ab}\cdot\mathbf{p}_b)}{m_a m_b}
          \right\}
         + \frac{1}{2}\sum_{a, b\neq c, c\neq a}\frac{m_a m_b m_c}{r_{ab}r_{ac}}.
```
