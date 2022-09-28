# The Euler Method

## A Simple Differential Equation

When learning numerical methods for solving equations, it's often
best to begin with simple problems with known solutions.  Then we
can easily see the strengths and weaknesses of different numerical methods.

We can begin with the simple differential equation
```math
\frac{dy}{dx} = \lambda y,
```
where ``\lambda`` is a constant.  The solution to this equation is
```math
y = e^{\lambda x},
```
which can be found directly by integration.


## A Discrete Solution

Numerical methods for differential equations are based on the idea
of discretization, where continuous functions become discrete functions
defined only at specific points.  Operators, like the derivative or
integral, are then approximated by algebraic operators, that can be
easily programmed for a computer.

To discretize this differential equation, we define the function
``y(x)`` at discrete points ``x_k``.  Here ``k`` is an index that
labels the points.  For simplicity, we will use a set of points with
uniform spacing, ``h``, defined as
```math
x_k = x_0 + k h
```
for some constant ``x_0``.

The derivative is formally defined as 
```math
\frac{dy}{dx} \equiv \lim_{\Delta x\to 0} \frac{y(x+\Delta x) - y(x)}{\Delta x},
```
where the limit is taken as ``\Delta x`` goes to zero.  We can approximate
this operation for discrete functions by dropping the limit and using a 
small, but finite, value for ``\Delta x``.  The derivative could
then be approximated as
```math
\frac{dy}{dx} \simeq \frac{y(x+h) - y(x)}{h},
```
where ``h`` is some small number. We might assume that the smaller ``h``
is, the better that this discrete derivative approximates the real
derivative.  This turns out to be true.

Using this approximation of the derivative, we can now try solving the
differential equation, ``\frac{dy}{dx} = \lambda y.``  We first define
some notation.  As ``y`` is defined only at points ``x_k``, 
we let ``y_k = y(x_k).``  Then, we can write the differential equation
```math
\frac{y_{k+1} - y_k}{h} = \lambda y_k.
```
We can rewrite this as
```math
y_{k+1} = y_k + h\lambda y_k = (1+h\lambda)y_k.
```
If we know the solution at ``x_k``, we can find the solution at the 
next point ``x_{k+1}`` using this equation.

## Example

The Jupyter notebook [EulerMethod.ipynb](https://github.com/davidneilsen/numex.jl/blob/main/examples/EulerMethod.ipynb) demonstrates the
Euler method to solve the ODE

## Going Further

(1) Try solving some different ODEs with the Euler method.  For example, solve
```math
\frac{dy}{dx} = -x y,
```
on the domain ``x\in[0,5]`` with the initial value ``y(0)=2``.  The exact solution has the form
```math
y = A e^{-x^2/2}.
```

(2) Try solving a system of equations
```math
\begin{aligned}
\frac{du}{dx} &= 3u + 2y,\\
\frac{dy}{dx} &= 4u + y,
\end{aligned}
```
on the domain ``0\le t\le 1`` with the initial conditions ``u(0) = 0`` and ``y(0) = 1.``  Try different values of ``h``, such as ``h = 0.2, 0.1, 0.05, 0.025``.

(3) How does the error in the Euler method depend on the discrete step size ``h``?  *Hint.* Expand the function ``y(x+h)`` in a Taylor series about ``x``
```math
y(x+h) = y(x) + \frac{dy}{dx} h + \frac{d^2y}{dx^2} h^2.
```


