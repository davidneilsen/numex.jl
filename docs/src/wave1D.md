# The Classical Wave Equation in 1D

The classical wave equation one dimension and Cartesian coordinates is
```math
\frac{\partial^2 \phi}{\partial^2 t} - c^2\frac{\partial^2 \phi}{\partial^2 x} = 0,
```
where ``c`` is the wave speed. To simplify the notation, 
we write define the notation that
```math
\phi_{,t} \equiv \frac{\partial \phi}{\partial t},
```
then the wave equation can be written as
```math
\phi_{,tt} - c^2\phi_{,xx} = 0.
```

## The Exact Solution

The one-dimensional wave equation has a simple analytic solution found by
d'Alembert.  If we introduce the new variables
```math
\xi = x - ct, \qqad \eta = x + ct,
```
the wave equation in these new variables simplifies to
```math
\phi_{,\xi\eta} = 0.
```
This equation has the general solution
```math
\phi = F(\xi) + G(\eta),
```
for arbitrary functions $F$ and $G$, or
```math
\phi(x,t) = F(x - ct) + G(x + ct).
```

## First-Order Variables for the Wave Equation

We first rewrite the second-order wave equation in first-order form.
Define the first order variables
```math
\begin{align}
\Pi &\equiv \phi_{,t}\\
\Phi &\equiv \phi_{,x},
\end{align}
```
then the wave equation can be written
```math
\Pi_{,t} = c^2\Phi_{,x}.
```
A second equation comes from the integration condition that 
mixed partial derivatives commute
```math
\frac{\partial^2 \phi}{\partial t \partial x} = \frac{\partial^2 \phi}{\partial x \partial t}.
```
In first-order form, this condition is
```math
\Phi_{,t} = \Pi_{,x}.
```
The system of equations that we will solve is
```math
\begin{align}
\Pi_{,t} &= c\Phi_{,x}\\
\Phi_{,t} &= \Pi_{,x}.
\end{align}
```

## Numerical Algorithm

There are several ways to solve the wave equation numerically.  We will
use the Method of Lines.  In this method, we first write the equations
as a semi-discrete system, where the spatial derivatives are discretized,
while the time derivatives remain analytic.  This gives a set of ordinary
differential equations (ODES) in time.

Let $\bf y$ be the 


