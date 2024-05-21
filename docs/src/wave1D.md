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

## Boundary Conditions

The wave equation is solved on a domain $a < x < b$.  Boundary conditions
must be given at the boundaries $x=a$ and $x=b$.  Sturm-Liouville theory
shows that solutions exist if the boundary conditions are 
Dirchlet, Neumann, Robin, or Periodic.

### Fixed Boundary Condition

The value of $\phi$ is specified at the boundary.  For example, a string
with a fixed endpoint corresponds to the boundary condition $phi(L,t) = 0$.

### Open Boundary Condition



### Sommerfeld Boundary Condition

Characteristic boundary condition.



## Numerical Algorithm

There are several ways to solve the wave equation numerically.  We will
use the Method of Lines.  In this method, we first write the equations
as a semi-discrete system, where the spatial derivatives are discretized,
while the time derivatives remain analytic.  This gives a set of ordinary
differential equations (ODES) in time.

Let $\bf y$ be the state vector
```math
{\bf y} = \begin{pmatrix} \Pi & \Phi \end{pmatrix}
```
and ${\bf f}({\bf y})$ be the RHS vector
```math
{\bf f}({\bf y}) = \begin{pmatrix} c\Phi & \Pi\end{pmatrix},
```
then the wave equation can be written
```math
\partial_t {\bf y} = \partial_x {\bf f}({\bf y}).
```



## Code

### Run in the Terminal

### Output

### Jupyter Notebook
