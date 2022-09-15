# Orbits in a Kerr Spacetime

The Kerr solution in general relativity describes spacetime geometry about 
a spinning point-mass of mass ``M``.  The angular momentum of the black hole
is given by the dimensionless angular momentum parameter ``a``, which has
values in the range ``0 \le a \le 1``.

The Kerr metric in Boyer-Lindquist coordinates is
```math
ds^2 = -\left(1-\frac{2r}{\Sigma}\right)\, dt^2
       - \frac{4ar\sin^2\theta}{\Sigma}\,dt\,d\phi
       + \sin^2\theta\left(r^2 + a^2 + \frac{2a^2 r \sin^2\theta}{\Sigma}\right)\,d\phi^2
       + \frac{\Sigma}{\Delta}\,dr^2
       + \Sigma\,d\theta^2,
```
where
```math
\begin{aligned}
\Sigma &\equiv r^2 + a^2\cos^2\theta,\\
\Delta &= r^2 - 2r + a^2.
\end{aligned}
```

Geodesics in the Kerr spacetime have three integrals of motion: the
orbital energy ``E``, the ``z``-component of the orbital angular
momentum ``L_z``, and the Carter constant ``Q``.  The geodesic equations
can then be written in terms of the particle's proper time ``\tau`` as
```math
\begin{aligned}
\frac{dr}{d\tau} &= \frac{1}{\Sigma} \pm \sqrt{R}\\
\frac{d\theta}{d\tau} &= \frac{1}{\Sigma} \pm \sqrt{\Theta}\\
\frac{d\phi}{d\tau} &= \frac{1}{\Sigma}\left[\frac{a}{\Delta}(2rE - aL_z) 
                       + \frac{L_z}{\sin^2\theta}\right]\\
\frac{dt}{d\tau} &= \frac{1}{\Sigma}\left[\frac{(r^2+a^2)^2E - 2arL_z}{\Delta}
                       - a^2 E \sin^2\theta \right].
\end{aligned}
```
Here ``R(r)`` and ``\Theta(\theta)`` are "quasi-potentials"  
```math
\begin{aligned}
R(r) &= -(1-E^2)r^4 + 2r^3 - \left[a^2(1-E^2) + L_z^2\right]r^2
        + 2(aE-L_z)^2 r - Q\Delta ,\\
\Theta(\theta) &= Q - \cos^2\theta\left\{a^2(1-E^2) + \frac{L_z^2}{\sin^2\theta}\right\}.
\end{aligned}
```

The equations can be simplified by defining the *Mino time* ``\lambda`` as
```math
\frac{d\tau}{d\lambda} \equiv \Sigma 
               \rightarrow 
\frac{d}{d\lambda} = \Sigma\frac{d}{d\tau}.
```

The geodesic equations of motion in a Kerr spacetime can be written in 
Mino time as
```math
\begin{aligned}
\dot r &= \Delta p_r\\
\dot p_r &= \frac{1}{2}\left\{ -\Delta' p_r^2 + \left(\frac{\Delta' R - \Delta R'}{\Delta^2}\right)\right\}\\
\dot \theta &= p_\theta\\
\dot p_\theta &= \frac{1}{2} \Theta^\theta\\
\dot \phi &= -\frac{\partial}{\partial L} \left(\frac{R}{\Delta} + \Theta\right)\\
\dot p_\phi &= 0\\
\dot t &= \frac{1}{2}\frac{\partial}{\partial E} \left(\frac{R}{\Delta} + \Theta\right)\\
\dot p_t &= 0.
\end{aligned}
```
Here the prime (``'``) and ``\theta`` superscripts indicate differentiation 
with respect to ``r`` and ``\theta``, respectively.

Reference: Gabe Perez-Giz, *From Measure Zero to Measure Hero: Periodic Kerr
Orbits and Gravitational Wave Physics,* (PhD Dissertation, Columbia University, 2011).
