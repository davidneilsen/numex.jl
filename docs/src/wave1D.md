# The Classical Wave Equation in 1D

The classical wave equation one dimension and Cartesian coordinates is
```math
\frac{\partial^2 \phi}{\partial^2 t} - c^2\frac{\partial^2 \phi}{\partial^2 x} = 0,
```
where ``c`` is the wave speed. To simplify the notation, 
we write define the notation that
```math
\phi_{,t} \equiv \frac{\partial \phi}{\partial t}.
```

To solve this equation numerically, we rewrite this 
second-order equation as first-order differential equations.  
Let's define the first order variables
```math
\begin{align}
\Pi &= \phi_{,t}\\
\Phi &= \phi_{,x},
\end{align}
```
then the wave equation can be written
```math
\Pi_{,t} = c^2\Phi_{,x}.
```
A second equation comes from the integration condition that mixed partial derivatives
commute
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


