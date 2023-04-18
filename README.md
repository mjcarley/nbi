NBI is a library and codes for the Nystrom Boundary Integral method.

It contains the tinyexpr parser and evaluation library written by
Lewis Van Winkle and released under the zlib licence. It is maintained
here:

https://github.com/codeplea/tinyexpr

* Prerequisites

You will need to have installed the following libraries:

https://github.com/mjcarley/blaswrap
https://github.com/mjcarley/sqt
https://github.com/mjcarley/wbfmm
https://github.com/mjcarley/agg

You should also install PETSc if you want to use it in place of the
basic built-in iterative solver:

https://petsc.org/

and make sure the environment variable PETSC_DIR is set and visible in
the shell.

* Installation

If you have downloaded the source from github or equivalent, you will
need the autotools suite and you generate the configure script with

. autogen.sh

To configure and install the code,

  ./configure [OPTIONS]
  make
  make install

For information on options to control configuration,

  ./configure --help
