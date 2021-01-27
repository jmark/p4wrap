
# This file is part of p4est
# Use `include /path/to/Makefile.p4est.mk' in your Makefile
# to use p4est in your project without autotools

p4est_prefix = /home/jmark/projects/p4wrap/sites/gcc/p4est/_prefix_
p4est_exec_prefix = ${p4est_prefix}
p4est_sysconfdir = ${p4est_prefix}/etc

include ${p4est_sysconfdir}/Makefile.sc.mk

# P4EST_CC and P4EST_CFLAGS may not be necessary for your project
P4EST_CC = mpicc
P4EST_CFLAGS = -g -O2

# These pull in p4est but none of its dependencies
P4EST_PKG_CPPFLAGS = -I${p4est_prefix}/include
P4EST_PKG_LDFLAGS = -L${p4est_exec_prefix}/lib
P4EST_PKG_LIBS = -lp4est

# These pull in everything needed by p4est
P4EST_CPPFLAGS =  $(SC_PKG_CPPFLAGS) $(P4EST_PKG_CPPFLAGS)
P4EST_LDFLAGS =  $(SC_PKG_LDFLAGS) $(P4EST_PKG_LDFLAGS)
P4EST_LIBS = $(P4EST_PKG_LIBS) $(SC_PKG_LIBS) -llapack -lblas -lz -lm   
