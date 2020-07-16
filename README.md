# P4WRAP v0.9

## Fortran 2008 wrapper around the 'p4est' [1] library.

** This project is still alpha and is subject to any changes at any time. **

Author: Johannes Markert <johannes.markert@jmark.de>

## Dependencies
* recent C99 and Fortran 2008 compilers
* installed 'MPI' library with compiler wrappers
* installs 'p4est' library automatically. See [1] for details.

## Compilation
```bash
  SHELL-PROMPT> make -C sites/[gcc|intel]
```

-------------------------------------------------------------------------------

[1] http://p4est.github.io
