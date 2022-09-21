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
=  \sum_\limits{j=1 \\ j\neq i}^{\infty} 
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

We can solve the ODEs

## Examples

