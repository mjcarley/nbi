#! /bin/sh

aclocal -I . \
&& automake --add-missing \
&& autoconf
