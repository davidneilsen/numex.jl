# White Dwarfs and Neutron Stars

Between 1930 and 1935, Subrahmanyan Chandrasekhar showed that white dwarf
stars have a maximum mass.  Stars beyond this mass limit become unstable
to collapse.  Chandrasekhar modeled the stars as degenerate Fermi gases.
In the non-relativistic limit, the pressure is proportional to the mass
density with an exponent of ``5/3``
```math
P = k \rho^{5/3}.\qquad {\rm (non-relativistic)}
```
This simple form for the equation of state is called a polytrope.
As the pressure in the star increases, however, the thermal motion of
the electrons reaches relativistic speeds.  In the relativistic limit, 
the equation of state softens to 
```math
P = k \rho^{4/3}.\qquad {\rm (relativistic)}
```
This equation of state has an upper mass limit of ``1.4~M_\odot`` for
the degenerate electron gas, now called the Chandrasekhar limit.
We now know that cold stars with masses that exceed this limit collapse
to neutron stars.

The fully relativistic solution for spherical stars was found by Oppenheimer
and Volkoff in 1939, building on the work of Tolman.  This solution also
has a limiting mass, and cold stars above this mass limit collapse to
form black holes.  The exact mass limit is difficult to determine, because
of our incomplete understanding of the nuclear equation of state at
high densities.  A very general analysis places an upper mass limit at 
``3~M_\odot``.  A combination of theoretical work and observational data
suggests that the limit could be as low as ``2.2~M_\odot.``  The mass
limit for highly spinning stars could be about 20% larger.

## The TOV Equations

The TOV equations are equations for a spherical star in
hydrostatic equilibrium in general relativity.  The pressure in the 
star is provided by the degeneracy pressure of an ideal Fermionic
gas.  This gas is characterized by the number density of 
particles ``n``, the mass density ``\rho``, the energy ``\epsilon``,
and the pressure ``P.``  Finally, ``m(r)`` gives the 
mass of the star contained within a radius ``r``, so that ``m(0) = 0,``
and ``\Phi`` a gravitational potential that specifies the spacetime
geometry.
The TOV equations are then
```math
\begin{aligned}
\frac{dm}{dr} &= 4\pi r^2 \epsilon\\
\frac{dP}{dr} &= - (\epsilon + P)\frac{m + 4\pi r^3 P}{r(r-2m)}\\
\frac{d\Phi}{dr} &= - \frac{1}{\epsilon + P}\frac{dP}{dr}
\end{aligned}
```

## Numerical Solutions

### Dimensionless Variables

### Example

# References

There are several good references on the TOV equations and 
the equations of state.  For now, these can be consulted for more
information.

1. Aaron Smith, "Tolman--Oppenheimer--Volkoff (TOV) Stars," (2012) [PDF](http://www.as.utexas.edu/astronomy/education/spring13/bromm/secure/TOC_Supplement.pdf).
2. Richard R. Silbar and Sanjay Reddy, "Neutron stars for undergraduates," *Am. J. Physics* **72** 892 (2004); [DOI](https://doi.org/10.1119/1.1703544), [arXiv](https://arxiv.org/abs/nucl-th/0309041).
3. Irina Sagert, Matthias Hempel, Carsten Greiner, Juergen Schaffner-Bielich, "Compact Stars for Undergraduates," *Eur. J. Phys.* **27** 577--610 (2006); [DOI](https://doi.org/10.1088/0143-0807/27/3/012), [arXiv](https://arxiv.org/abs/astro-ph/0506417).
