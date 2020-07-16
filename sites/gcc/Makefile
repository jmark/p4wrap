SOURCEDIR = ../../source
BUILDDIR = build

P4EST_PREFIX="$(shell pwd)/p4est/_prefix_"

CC = mpicc

CCFLAGS =

CCFLAGS += -O3
#CCFLAGS += -g

CCFLAGS += -fPIC
CCFLAGS += -std=c99
CCFLAGS += -pedantic

CCFLAGS += -Wall
CCFLAGS += -Werror
CCFLAGS += -Wno-array-bounds
CCFLAGS += -Wno-unused-variable
CCFLAGS += -Wno-unused-but-set-variable
CCFLAGS += -Wno-unused-function

CCFLAGS += -I$(SOURCEDIR)
CCFLAGS += -Ip4est/_prefix_/include

.PHONY: default wrappers staticlib clean clear

default: staticlib

clean clear:
	rm -rf p4est/*
	rm -rf $(BUILDDIR)

p4est/.dummy:
	git submodule update --init --recursive
	cd p4est && ./bootstrap
	cd p4est && ./configure --enable-mpi --enable-debug --prefix="$(P4EST_PREFIX)"
	cd p4est && make -j4
	cd p4est && make install
	touch $@

wrappers: $(BUILDDIR)/p4wrap.o $(BUILDDIR)/p8wrap.o

$(BUILDDIR)/p4wrap.o: $(SOURCEDIR)/p4wrap.c $(SOURCEDIR)/p4wrap.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c -o $@ $<

$(BUILDDIR)/p8wrap.o: $(SOURCEDIR)/p8wrap.c $(SOURCEDIR)/p4wrap.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c -o $@ $<

staticlib: $(BUILDDIR)/libp4wrap.a

$(BUILDDIR)/libp4wrap.a: p4est/.dummy wrappers
	ar rcs $@ $(shell find p4est -type f -iname '*.o') $(shell find build -type f -iname '*.o')
