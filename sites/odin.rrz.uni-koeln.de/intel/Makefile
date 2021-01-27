SOURCEDIR = ../../../source
BUILDDIR = build

CC = mpiicc

CCFLAGS =
CCFLAGS += -O3 
CCFLAGS += -fPIC
CCFLAGS += -std=c99

CCFLAGS += -I$(SOURCEDIR)

.PHONY: default wrappers staticlib clean clear

default: staticlib

clean clear:
	rm -rf $(BUILDDIR)

staticlib: $(BUILDDIR)/libp4wrap.a
wrappers: $(BUILDDIR)/p4wrap.o $(BUILDDIR)/p8wrap.o

$(BUILDDIR)/p4wrap.o: $(SOURCEDIR)/p4wrap.c $(SOURCEDIR)/p4wrap.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c -o $@ $<

$(BUILDDIR)/p8wrap.o: $(SOURCEDIR)/p8wrap.c $(SOURCEDIR)/p4wrap.h
	mkdir -p build
	$(CC) $(CCFLAGS) -c -o $@ $<

$(BUILDDIR)/libp4wrap.a: wrappers
	ar rcs $@ $(shell find build -type f -iname '*.o')
