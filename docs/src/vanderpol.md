# The van der Pol Oscillator

The van der Pol oscillator is
```math
\frac{d^2x}{dt^2} - \mu(1-x^2)\frac{dx}{dt} + x = 0.
```
Here ``\mu`` is a scalar parameter that controls the strength of the 
nonlinear damping.

Writing this in first-order form
```math
\begin{aligned}
\dot x &= y\\
\dot y &= \mu(1-x^2)y - x
\end{aligned}
```

## The Hamiltonian

The canonical momentae  are 
```math
\begin{aligned}
p_x &= \dot y + \mu(1-x^2)y\\
p_y &= \dot x.
\end{aligned}
```
The Hamiltonian is
```math
H = p_x p_y + xy -\mu(1-x^2)y p_y.
```
Hamilton's equations of motion are
