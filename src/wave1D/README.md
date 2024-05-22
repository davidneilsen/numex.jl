# The Classical Wave Equation in 1D


Date: *Wed May 22 09:21:19 MDT 2024*

This code was able to use VTK output in the past.  Apparently, the VTKWrite.jl package has changed, so that it no longer supports 1D data.  I changed the output to write the curve data format used by VisIT.  VisIT can be used to visualize the time-dependent data.


Date: *Fri Sep  9 13:08:02 MDT 2022*

Begin with the Jupyter notebook `WaveEqSolver.ipynb`.  The result of the simulation is saved as `Pi.gif`.

The wave equation solver can also be run directly from the commandline as

    julia solver.jl <nsteps>

The animated output files are `aPi.gif` and `aPhi.gif`.

