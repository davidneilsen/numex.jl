# Maxwell Equations

*Finite Difference Time Domain Maxwell Equations*

The simple example solves the Maxwell equations in two dimensions.

## Equations of motion

The Maxwell equations are
```math
\nabla \times E = -\frac{\partial B}{\partial t}.
```

```math
\begin{aligned}
\nabla\cdot\mathbf{E}  &= 4 \pi \rho \\
\nabla\cdot\mathbf{B}  &= 0 \\
\nabla\times\mathbf{E} &= - \frac{1}{c} \frac{\partial\mathbf{B}}{\partial t} \\
\nabla\times\mathbf{B} &= - \frac{1}{c} \left(4 \pi \mathbf{J} + \frac{\partial\mathbf{E}}{\partial t} \right)
\end{aligned}
```

## Numerical method

## Try it out


