* Examples for use of NBI in Laplace and Helmholtz problems

The examples in each directory are run using the scripts
solve-helmholtz and solve-laplace, which can also be used to configure
the problem. Each directory contains scripts make-geometry and
make-visualization-grid which generate the surface mesh and a
visualization mesh, respectively. There are also configurable boundary
condition (source) files which are used in running the tests.

The test cases are intended to be run with a source contained inside
the surface. When the boundary condition is generated in this way, the
resulting solution should be identical to the potential computed using
the interior source. The NBI solver returns error estimates (L_2,
L_inf, and r.m.s.) computed as the difference between the potential
specified in the boundary condition file, and that computed by the
solver.

The output from the calculation which be visualized using GMSH is two
.msh files, one named according to the geometry of the problem,
e.g. sphere.msh, and one called grid.msh which shows the computed
field around the surface. 

To run all of the test cases, with the 

SphereLaplace: solve for the Laplace potential on a sphere under
	       excitation from an internal point source. To run the
	       example in the directory,

	       ./sphere-laplace

SphereHelmholtz: solve for the Helmholtz (wave equation) potential on
		 a sphere under excitation from a sinusoidally varying
		 ring source (this can be used to simulate rotor and
		 propeller noise). To run the example in the
		 directory,

		 ./sphere-helmholtz

Scattering: solve for scattering of a plane wave by a sphere, with
	    comparison to the analytical series solution. To run the
	    example in the directory,

	    ./sphere-scattering

	    The README file in the directory contains instructions on
	    comparing the scattered field to an analytical calculation
	    using Octave. 

Rotor: solve the problem of scattering of noise from an azimuthally
       varying circular source, which simulates a rotor interacting
       with a simple nacelle

Catseye: solve for a point source inside a cat's eye geometry,
	 generated using a GMSH input file.  To run the example in the
	 directory,

	 ./catseye-laplace

* Parameters

The shell scripts take a number of options, which can be set from the
command line using the following options:

-k k:        acoustic wavenumber, for Helmholtz problems
-a nqa:      number of quadrature points in rule used for adaptive integration
-p nqp:      number of quadrature points on discretized patches
-u nqu:      number of quadrature points on upsampled patches
-t tol:      adaptive integration tolerance
-o order:    order of singular integration
-e eta:      near/far field cutoff parameter
-d depth:    maximum recursion depth for adaptive quadrature
-T nthreads: number of threads to use in threaded computations
-s src:      boundary condition (source) filename