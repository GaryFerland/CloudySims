.PHONY: all check64 check32 clean

CXX = g++
CXXFLAGS = -g -W -Wall -ansi -std=c++11 -O0 -I../src -I../cpuid -I../unittest-cpp
LDFLAGS = -lUnitTest++ -lvhsum
SED = sed
MV = mv
CCDEP = $(CXX)
CCDEPFLAGS = $(CXXFLAGS)

VHSUMLOC32 = $(shell realpath ../lib32)
VHSUMLOC64 = $(shell realpath ../lib64)

ifeq ($(MAKECMDGOALS),check32)
  CXXFLAGS += -m32
  LDFLAGS += -L../unittest-cpp/lib32 -L../lib32 -m32 -Wl,-rpath,$(VHSUMLOC32)
  SHAREDLIB = ../lib32/libvhsum.so
else
  LDFLAGS += -L../unittest-cpp/lib64 -L../lib64 -Wl,-rpath,$(VHSUMLOC64)
  SHAREDLIB = ../lib64/libvhsum.so
endif

test_src = TestMain.cc TestCore.cc TestScalar.cc TestSSE2.cc TestAVX2.cc TestAVX512f.cc
test_obj = $(patsubst %.cc, %.o, $(test_src))
test_deps = $(patsubst %.cc, %.d, $(test_src))

define make-depend
  $(CCDEP) $(CCDEPFLAGS) -MM $1 | \
  $(SED) -e 's,\($(notdir $2)\) *:,$(dir $2)\1: ,' > $3.tmp
  $(SED) -e 's/#.*//' \
      -e 's/^[^:]*: *//' \
      -e 's/ *\\$$//' \
      -e '/^$$/ d' \
      -e 's/$$/ :/' $3.tmp >> $3.tmp
  $(MV) $3.tmp $3
endef


all: check64

check64: TestMain run_tests.sh
	@echo "=============================="
	@echo " Starting VectorHash tests... "
	@echo "=============================="
	./TestMain
	./run_tests.sh

check32: TestMain
	@echo "=============================="
	@echo " Starting VectorHash tests... "
	@echo "=============================="
	./TestMain

clean:
	rm -f *.o
	rm -f *.d
	rm -f TestMain

%.o: %.cc
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CXXFLAGS) -c -o $@ $<

TestMain: $(test_obj) $(SHAREDLIB)
	$(CXX) $(test_obj) $(LDFLAGS) -o TestMain

# This Makefile should only be called from the Makefile in the parent directory
# That Makefile will have built these libs already, so nothing needs to be done here...
../lib32/libvhsum.so:

../lib64/libvhsum.so:

ifneq "$(MAKECMDGOALS)" "clean"
-include $(test_deps)
endif
