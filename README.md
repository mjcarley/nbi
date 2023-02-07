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

There appears to be an optimization bug in gcc version 11 (at least),
which seriously affects the performance of the code, so it is
recommended to switch off optimization of NBI (most of the calculation
is done using other libraries, so this does not affect computation
times too badly), using

  CFLAGS="-O0 -g" ./configure ...

