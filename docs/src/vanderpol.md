# The van der Pol Oscillator

The van der Pol oscillator is a nonlinear oscillator whose solutions
can show deterministic chaos.  The oscillator was originally derived
to study circuits with vacuum tubes, but it has found applications
in biophysics, geology, and in modeling vocal folds in acoustics.

The oscillator satisfies the equation
```math
\frac{d^2x}{dt^2} - \mu(1-x^2)\frac{dx}{dt} + x = 0.
```
The middle term proportional to the velocity ``dx/dt`` is a 
damping term that depends on the scalar parameter ``\mu``.  
The dependence on ``x^2`` makes this damping nonlinear.

While the linear damped oscillator can be solved analytically,
the nonlinear oscillator does not have a known analytic solution.
To solve the system numerically, we first write it in first-order
form by introducing the new variable ``y = dx/dt``
```math
\begin{aligned}
\dot x &= y\\
\dot y &= \mu(1-x^2)y - x
\end{aligned}
```
Again, we are using an overdot (``\dot{}``) to represent a time
derivative ``d/dt``.  These equations can now be solved
by many different ODE solvers.  

An example of solving the van der Pol equations is found in this notebook
[vanderPol.ipynb](https://github.com/davidneilsen/numex.jl/blob/main/examples/vanderPol.ipynb).  In this example, I use the `DifferentialEquations.jl` 
package to solve the equations.

# Explore the Limit Cycle

One method for studying dynamical systems is to plot the solutions in phase space.
A phase space plot for a simple system is made by plotting the the position and
momentum at regular intervals on a two-dimensional plot.  For example, the phase
space plot for a simple harmonic oscillator would simple be a closed loop.

In some nonlinear systems, there is an attractor or a closed orbit in phase space
towards which the system tends to evolve. If the system begins in a state away
from the attractor, it tends to evolve towards the attracting solution.  The 
van der Pol system has an attractor.  The example file [vanderPol.ipynb](https://github.com/davidneilsen/numex.jl/blob/main/examples/vanderPol.ipynb) solves the ODEs for 100
different, randomly chosen sets of initial conditions and then plots the velocity
``\dot x`` vs. the position ``x``.  We see that most of the solutions evolve to the
attractor solution, while a few simple stop oscillating (the dots on the ``x``-axis
with ``\dot x = 0``).
