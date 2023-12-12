#! /bin/sh

mkdir -p m4

aclocal -I .
autoheader
automake --add-missing
autoconf
