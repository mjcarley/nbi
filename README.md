NBI is a library and codes for the Nystrom Boundary Integral method.

It contains the tinyexpr parser and evaluation library written by
Lewis Van Winkle and released under the zlib licence. It is maintained
here:

https://github.com/codeplea/tinyexpr

# Prerequisites

You will need to have installed the following libraries:

- https://github.com/mjcarley/blaswrap
- https://github.com/mjcarley/sqt
- https://github.com/mjcarley/wbfmm
- https://github.com/mjcarley/agg

If you clone the github repository using

`git clone --recursive https://github.com/mjcarley/nbi`

these libraries will be downloaded automatically in the `externals`
subdirectory of `nbi`. 

You should also install PETSc if you want to use it in place of the
(very) basic built-in iterative solver:

https://petsc.org/

and make sure the environment variable PETSC_DIR is set and visible in
the shell.

# Installation

If you have downloaded the source from github or equivalent, you will
need the autotools suite. Generate the configure script with

`. autogen.sh`

To configure and install the code,

`./configure [OPTIONS]`

`make`

`make install`

For information on options to control configuration, including the
installation location and where to find any required libraries:

  `./configure --help`

# Using and testing NBI

The main documentation, including details of examples, is in the
`.../doc` subdirectory of the distribution. The `.../examples`
subdirectory contains a number of examples which can be run using the
scripts test-all-laplace and test-all-helmholtz. Further details are
given in `.../examples/README.md`
