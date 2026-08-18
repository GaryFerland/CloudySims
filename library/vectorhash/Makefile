.PHONY: all default lib32 lib64 testclean clean distclean check check32 install install-lib

# on modern Mac Xcode and FreeBSD releases, g++ is disappearing in favor of clang++
# so make sure we have a fallback in case g++ is missing...
DEFCXX := ${shell which g++ 2> /dev/null}
ifeq ($(DEFCXX),)
  DEFCXX := ${shell which clang++ 2> /dev/null}
endif

override CXX = $(DEFCXX)
override CXXFLAGS = -g -W -Wall -Wno-unused-command-line-argument -ansi -std=c++11 -O3 -funroll-loops -fPIC -fvisibility=hidden

INSTALLDIR = /usr/local

GZIP := ${shell which gzip}
ifneq ($(GZIP),)
  GZIP = gzip
  GZEXT = .gz
else
  GZIP = touch
  GZEXT =
endif

# DEP_GOALS will be an empty string for targets that don't need Makefile.dep
NODEP_GOALS = testclean clean distclean install
DEP_GOALS = $(filter-out $(NODEP_GOALS),$(MAKECMDGOALS))

# the default target does need Makefile.dep
ifeq ($(MAKECMDGOALS),)
  DEP_GOALS = need_dependencies
endif

ifneq ($(DEP_GOALS),)
  $(shell test ! -f Makefile.dep && script/generate_deps.sh "$(CXX)" "$(CXXFLAGS)" > Makefile.dep)
endif

all: default

-include Makefile.dep

default: bin/vh128sum bin/vh256sum bin/vh512sum $(LIB64NAME)

testclean:
	cd tests; \
	$(MAKE) clean

clean: testclean
	rm -f Makefile.dep
	rm -f lib64/*.o
	rm -f lib64/*/*.o
	rm -f lib32/*.o
	rm -f lib32/*/*.o
	rm -f lib64/libvhsum.*
	rm -f lib32/libvhsum.*
	rm -f bin/vh128sum
	rm -f bin/vh256sum
	rm -f bin/vh512sum

distclean: clean
	cd unittest-cpp; \
	$(MAKE) clean

bin/vh128sum: lib64/vectorhash.o $(LIB64NAME)
	$(CXX) lib64/vectorhash.o $(LDFLAGS) -o $@

bin/vh256sum: bin/vh128sum
	ln -f bin/vh128sum bin/vh256sum

bin/vh512sum: bin/vh128sum
	ln -f bin/vh128sum bin/vh512sum

unittest-cpp/lib64/libUnitTest++.a:
	cd unittest-cpp; \
	$(MAKE)

unittest-cpp/lib32/libUnitTest++.a:
	cd unittest-cpp; \
	$(MAKE) lib32

check: unittest-cpp/lib64/libUnitTest++.a $(LIB64NAME)
	cd tests; \
	$(MAKE)

check32: unittest-cpp/lib32/libUnitTest++.a $(LIB32NAME)
	cd tests; \
	$(MAKE) check32

install: install-lib
	mkdir -p $(INSTALLDIR)/bin
	cp -af bin/vh*sum $(INSTALLDIR)/bin
	mkdir -p $(INSTALLDIR)/include
	cp -af src/vectorhash.h $(INSTALLDIR)/include
	mkdir -p $(INSTALLDIR)/man/man1
	cp -af man/vh128sum.1 $(INSTALLDIR)/man/man1
	cd $(INSTALLDIR)/man/man1; \
	$(GZIP) -f vh128sum.1; \
	ln -sf vh128sum.1$(GZEXT) vh256sum.1$(GZEXT); \
	ln -sf vh128sum.1$(GZEXT) vh512sum.1$(GZEXT)
	mkdir -p $(INSTALLDIR)/man/man3
	cp -af man/VectorHash.3 $(INSTALLDIR)/man/man3
	cd $(INSTALLDIR)/man/man3; \
	$(GZIP) -f VectorHash.3; \
	ln -sf VectorHash.3$(GZEXT) VectorHashSIMD.3$(GZEXT)
