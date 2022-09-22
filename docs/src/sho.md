# The Simple Harmonic Oscillator


As mentioned before, it is good to test numerical methods
with problems with known solutions.  In this example we will use
the simple harmonic oscillator.

We will build on the previous example in two ways.  First, we 
solve a system of differential equations, and we will improve
on Euler method.

## The Differential Equation and Solution

Let's now consider the
simple harmonic oscillator
```math
\frac{d^2 y}{dt^2} = -\omega^2 y.
```
The equation is solved on a domain $0 \le t$ as an initial value problem,
with ``y(0)`` and ``\dot{y}(0)`` specified at the initial time.
The solution of this equation can be written
```math
y(t) = A \sin(\omega t) + B\cos(\omega t),
```
where the constants `A` and `B` are chosen to match the
initial conditions.  For simplicity, we use the notation that a dot
over a variable represents a time derivative
```math
\dot y \equiv \frac{dy}{dt}.
```

## Writing as a First-Order System

The simple harmonic oscillator equation (SHO) is a second-order equation.
Many numerical methods for the initial value problem have been
developed for first-order equations.  We can rewrite the SHO as a 
system of first-order equations by introducing a new variable ``u = \dot y``.
Then
```math
\begin{aligned}
\dot y &= u\\
\dot u &= -\omega^2 y.
\end{aligned}
```


## The Midpoint Method

While the Euler method is very simple to use, the method is not very accurate.
The Midpoint method is a similar method, but has higher accuracy.  
Consider the differential equation
```math
\frac{dy}{dt} = f(t,y).
```
The Euler method approximates the derivative ``dy/dt`` as a constant ``f(t_i)``
over the domain ``[t_i,t_{i+1}]``, or
```math
y_{i+1} = y_i   + \Delta t f(t_i,y_i)
```
A better approximation might be to use the value of the derivative at the
midpoint of ``t_i`` and ``t_{i+1}``, rather than at the endpoint ``t_i``,
so that
```math
y_{i+1} = y_i   + \Delta t f(t_{i+1/2}, y_{i+1/2})
```
But how do we find the value of ``y_{i+1/2}``?  We could use an Euler
step of ``\Delta t/2``.  Thus, the midpoint method is a two-step method.
We first approximate the solution at ``t_{i+1/2}`` using the Euler method,
then we go back to ``t_i`` and take a full step with the derivative
evaluated at ``t_{i+1/2}``
```math
\begin{aligned}
y_{i+1/2} &=  y_i   + \frac{1}{2}\Delta t f(t_i,y_i)\\
y_{i+1} &=  y_i   + \frac{1}{2}\Delta t f(t_{i+1/2},y_{i+1/2}).
\end{aligned}
```

## Other Second-Order Methods

The Midpoint method is not the only integrator with second-order accuracy.
There is a larger class of Runge-Kutta integrators of different orders,
and this is just one of the second-order Runge-Kutta methods, sometimes 
called RK2 methods.  A second RK2 method is
```math
\begin{aligned}
y^{\star} &=  y_i   + \Delta t f(t_i,y_i)\\
y_{i+1} &=  y_i   + \frac{1}{2}\Delta t 
\left( f(t_i,y_i) + f(t_{i+1},y^{\star})\right).
\end{aligned}
```
Beginning with the derivative at the left end-point, ``f(t_i, y_i)``,
the Euler method is used to approximate the derivative at the right
end-point, ``f(t_{i+1},y^{\star})``.  We see that ``y^\star`` is
a first approximation to ``y_{i+1}``, and this step is called the
*predictor* step.
The approximation to ``y_{i+1}`` is refined by stepping from
``t_i`` to ``t_{i+1}`` using the average of the two derivatives.
This is the *corrector* step of the algorithm.  


## Convergence

The Euler method is only first-order accurate, which means that the 
error in the
discrete solution is proportional to the step-size ``\Delta t``.
This is often written as ``{\rm O}(\Delta t)``.
If the step size is reduced by a factor of 2, then the error also 
drops by factor of 2.  However, the higher-resolution solution requires
even more work to compute, as it doubles the number of steps to reach
the final time.  

The Midpoint method is second-order accurate, or ``{\rm O}(\Delta t^2``.
If the step size is reduced by a factor or 2, then the error drops by
a factor of ``2^2=4``.  As the higher-order solution typically only
doubles the amount of work done, this is advantageous because the error 
drops faster than the increase in workload.  
Even higher order methods could be advantageous by this analysis. 
The fourth-order Runge-Kutta method, with ``{\rm O}(\Delta t^2)``,
is widely used for many conventional ODE systems.

## The DifferentialEquations.jl Package

There are many different numerical algorithms for solving ODEs.  Some
are very simple, such as the Euler and RK2 methods seen here.  There
are many types of differential equations, however, and many require
specialized numerical solvers.  These solvers can also be very complicated.
Thus, rather than writing our own solvers for each problem, it is good
to find reliable packages of ODE solvers.

Julia has one of the best packages for solving ordinary differential 
equations called [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). 


## Example

Consider a driven oscillator with the equations of motion written
in first-order form 
```math
\begin{aligned}
\dot y &= u\\
\dot u &= -100y + 8\cos(9t).
\end{aligned}
```
This oscillator is driven at a frequency just under the natural frequency
``\omega=10``, so the solution will have beats.  The exact solution of
this system has the form
```math
y(t) = A \cos(10 t) + B \sin(10 t) + \frac{8}{19}\cos(9 t).
```
where the constants ``A`` and ``B`` are determined by initial conditions.


The Jupyter notebook [EulerMethod.ipynb](https://github.com/davidneilsen/numex.jl/blob/main/examples/EulerMethod.ipynb) demonstrates the
Euler method to solve the ODE

## Going Further


