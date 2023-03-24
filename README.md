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

You should also install PETSc if you want to solve any scattering
problems:

https://petsc.org/

and make sure the appropriate environment variables are set.

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
