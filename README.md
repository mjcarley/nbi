Autotools Template
==================

Overview
--------

This is a skeleton project for a source tree based on the [GNU Build System](http://www.gnu.org/software/automake/manual/html_node/GNU-Build-System.html). The steps below walk through the process of running autotools so that the project can build built using the standard way:

    ./configure && make && make install

Prerequisites
-------------

You will need to install the following GNU tools:

    autoconf
    automake


Step-by-Step
------------

Clone this repository:

    git clone https://github.com/gizero/autotools-skeleton.git

Generate the configure script:

    autoreconf -ivf

Configure and build the project:

    ./configure
    make
