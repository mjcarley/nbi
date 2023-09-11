---
title: 'NBI: A library for Nystrom Boundary Integral calculations'
tags:
  - Nystrom Boundary Integral method
  - Helmholtz equation
  - Laplace equation
  - Fast Multipole Method
  - acoustics
  - potential
authors:
  - name: Michael J. Carley
    orcid: 0000-0003-2965-8984
    equal-contrib: true
    affiliation: 1
	corresponding: true
affiliations:
 - name: Department of Mechanical Engineering, University of Bath, UK
   index: 1
date: 7 September 2023
bibliography: paper.bib
--- 

# Summary

Boundary integral methods are a popular approach to solving potential
problems, in particular for the Laplace and Helmholtz problems in such
areas as fluid dynamics, acoustics, and electromagnetism. They are a
natural choice for the solution of problems in unbounded domains, such
as wave scattering, where the radiation boundary condition is
automatically satisfied by the nature of the formulation. The Nystrom
method is one approach to the solution of boundary integral problems,
which lends itself to acceleration using the Fast Multipole Method
(FMM). NBI is a library for the solution of boundary integral problems
based on the approach of [@greengard-oneil-rachh-vico21] and the FMM
methods of [@gumerov-duraiswami03; @gumerov-duraiswami05;
@gumerov-duraiswami09]. The code includes a number of executables
which can be used to set up and solve problems on realistic
geometries, with a number of examples provided for testing of the
solver. Results can be visualized using GMSH [@geuzaine-remacle09], a
standard free meshing program.

# Statement of need

A number of free boundary integral codes exist which can be used for
problems of the type handled by NBI [@kirkup07; @betcke-scroggs21, for
example], with different approaches to discretization and solution,
ranging from collocation methods with direct solvers, suitable for
relatively small problems, to FMM-accelerated Galerkin techniques with
iterative solvers, which can be used on large complex geometries. NBI
implements recent work on high-order solution of boundary integral
problems [@greengard-oneil-rachh-vico21] in a form which can be used
for a range of engineering applications. The principal motivation is
the solution of potential problems in acoustics and aerodynamics. As
well as the library proper, there are codes for discretization of
general geometries supplied in a number of formats; matrix assembly
and problem solution; field evaluation; and post-processing and
visualization of results. A built-in parser gives a flexible and
intuitive means of evaluating boundary conditions, which are supplied
as analytical expressions parsed and evaluated on the surface. This
allows boundary conditions to be specified in a form which makes it
easier to perform parametric studies by systematically modifying
internal variables. 

# Mathematics

The Laplace and Helmholtz potentials generated by a surface source
distribution are given by a boundary integral
\begin{equation}
\phi(\mathbf{x})
	=
	\int_{S}
	\phi(\mathbf{y})\frac{\partial}{\partial n}G(\mathbf{x},\mathbf{y})
	-
	\frac{\partial\phi(\mathbf{y})}{\partial n}G(\mathbf{x},\mathbf{y})
	\,
	\mathrm{d}S(\mathbf{y}),
\end{equation}
where $G$ is the Green's function for the problem. For points lying on
the surface, the left hand side is replaced by $\phi/2$ and the
equation is interpreted as a boundary integral equation to be solved
subject to 

# Features

Details of the problem setup and solution procedure are given in the
code documentation, but the basic steps are:

- generation of a surface discretization;
- assembly of the system matrices;
- solution of the problem subject to a specified boundary condition;
- postprocessing, including evaluation of the potential field, and
  visualization.

## Surface generation and representation

## Quadratures

[@bremer-gimbutas13]

## Solving problems

## Postprocessing and visualization

![GMSH visualisation of output from sphere scattering
example.\label{fig:scattering}](scattering.png){width=50%}

Figure \autoref{fig:scattering} shows the GMSH visualisation of the
scattered field from a sphere subject to plane wave excitation, one of
the test cases included in the package. 

![GMSH visualisation of output from rotor noise
example.\label{rotor}](rotor.png){width=50%}

Figure \autoref{rotor} shows the GMSH visualisation of the
field scattered from an ellipsoid, subject to an incident field from a
ring source, a simple model for scattering of rotor noise by a
nacelle. 

# References