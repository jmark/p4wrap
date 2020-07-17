# P4WRAP v0.99

## Fortran 2008 wrapper around the 'p4est' [1] library.

** This project is still alpha and is subject to any changes at any time. **

Author: Johannes Markert <johannes.markert@jmark.de>

## Dependencies
* recent C99 and Fortran 2008 compilers
* installed 'MPI' library with compiler wrappers

## Compilation
```bash
# fetch 'p4est'
git submodule update --init --recursive

# run make (adapt Makefile to your system)
make -C sites/[gcc|intel]

# copy static library into your project
cp sites/[gcc|intel]/build/libp4wrap.a path/to/your/project/build/...

# copy Fortran wrappers into your project
cp wrapper/*.f90 path/to/your/project/source/...
```

Have fun!

-------------------------------------------------------------------------------

[1] http://p4est.github.io
