# The Classical Wave Equation in 1D

The classical wave equation one dimension and Cartesian coordinates is
```math
\frac{\partial^2 \phi}{\partial^2 t} - c^2\frac{\partial^2 \phi}{\partial^2 x} = 0,
```
where ``c`` is the wave speed. To simplify the notation, 
we write define the notation that
```math
\partial_t \phi \equiv \frac{\partial \phi}{\partial t},
```
then the wave equation can be written as
```math
\partial_{tt}\phi - c^2\partial_{xx}\phi = 0.
```

## The Exact Solution

The one-dimensional wave equation has a simple analytic solution found by
d'Alembert.  If we introduce the new variables
```math
\xi = x - ct, \qquad \eta = x + ct,
```
the wave equation in these new variables simplifies to
```math
\partial_{\xi}\partial_{\eta} \phi= 0.
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
\Pi &\equiv \partial_t\phi\\
\Phi &\equiv \partial_x\phi,
\end{align}
```
then the wave equation can be written
```math
\partial_t\Pi = c^2\partial_x\Phi.
```
A second equation comes from the integration condition that 
mixed partial derivatives commute
```math
\frac{\partial^2 \phi}{\partial t \partial x} = \frac{\partial^2 \phi}{\partial x \partial t}.
```
In first-order form, this condition is
```math
\partial_t\Phi = \partial_x\Pi.
```
The system of equations that we will solve is
```math
\begin{align}
\partial_t\Pi &= c\partial_x \Phi\\
\partial_t\Phi &= \partial_x \Pi.
\end{align}
```
Let $\bf y$ be the state vector
```math
{\bf y} = \begin{pmatrix} \Pi \\ \Phi \end{pmatrix}
```
and ${\bf f}({\bf y})$ be the RHS vector
```math
{\bf f}({\bf y}) = \begin{pmatrix} c\Phi \\ \Pi\end{pmatrix},
```
then the wave equation can be written
```math
\partial_t {\bf y} = \partial_x {\bf f}({\bf y}).
```


## Boundary Conditions

The wave equation is solved on a finite domain $a < x < b$, and 
boundary conditions (BCs) must be provided at the boundaries $x=a$ and $x=b$
to give a well-defined problem.  In general, we can provide boundary
information in four ways:  we can specify 
  1. the the value of $\phi$ on the boundary (Dirichlet BC), 
  2. the spatial derivative $\partial_x \phi$ (Neumann BC), 
  3. a linear combination of $\phi$ and $\partial_x \phi$ (Robin BC), or
  4. periodic boundary conditions.

### Fixed Boundary Condition

Assume that $\phi(x,t)$ represents the displacement of a string.  If 
end of the string is fixed at $x=a$, then the boundary condition is 
$\phi(a,t) = 0$.  This could be generalized to make the displacement of the
string a known function of time, $\phi(a,t) = \phi_L(t)$.


### Open Boundary Condition

If one end of the string is free or open, then the boundary condition
is $\partial_x \phi = 0$.

### Sommerfeld Boundary Condition

The Sommerfeld boundary condition is an out-going condition to allow the
wave to pass unchanged through the boundary.
The coordinates $\xi$ and $\eta$ are the right-moving and left-moving
characteristic directions.
We can impose the condition $\partial_\xi \phi = 0$ on the left, and 
$\partial_\eta \phi= 0$ on the right.


## Numerical Algorithm

There are several ways to solve the wave equation numerically.  We will
use the Method of Lines.  In this method, we first write the equations
as a semi-discrete system, where the spatial derivatives are discretized,
while the time derivatives remain analytic.  This gives a set of ordinary
differential equations (ODES) in time.



## Code

### Run in the Terminal

### Output

### Jupyter Notebook
