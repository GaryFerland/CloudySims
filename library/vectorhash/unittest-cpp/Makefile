CXX = g++
CCFLAGS = -std=c++11 -g -W
LDFLAGS =
LDLIBS = -L.
SED = sed
MV = mv
RM = rm
CCDEP = $(CXX)
CCDEPFLAGS = -std=c++11 -g -W

.PHONY: all lib32 clean

.SUFFIXES: .o .cpp

lib = libUnitTest++.a
lib64 = lib64/$(lib)
lib32 = lib32/$(lib)
test = TestUnitTest++

RANLIB := ${shell which ranlib}
ifneq ($(RANLIB),)
  RANLIB = ranlib
else
  RANLIB = ar ts
endif

src = UnitTest++/AssertException.cpp \
	UnitTest++/Checks.cpp \
	UnitTest++/CompositeTestReporter.cpp \
	UnitTest++/CurrentTest.cpp \
	UnitTest++/DeferredTestReporter.cpp \
	UnitTest++/DeferredTestResult.cpp \
	UnitTest++/MemoryOutStream.cpp \
	UnitTest++/ReportAssert.cpp \
	UnitTest++/RequiredCheckException.cpp \
	UnitTest++/RequiredCheckTestReporter.cpp \
	UnitTest++/Test.cpp \
	UnitTest++/TestDetails.cpp \
	UnitTest++/TestList.cpp \
	UnitTest++/TestReporter.cpp \
	UnitTest++/TestReporterStdout.cpp \
	UnitTest++/TestResults.cpp \
	UnitTest++/TestRunner.cpp \
	UnitTest++/ThrowingTestReporter.cpp \
	UnitTest++/TimeConstraint.cpp \
	UnitTest++/XmlTestReporter.cpp

ifeq ($(MSYSTEM), MINGW32)
	src += UnitTest++/Win32/TimeHelpers.cpp
else
	src += UnitTest++/Posix/SignalTranslator.cpp \
	UnitTest++/Posix/TimeHelpers.cpp
endif

test_src = UnitTest++/tests/Main.cpp \
	UnitTest++/tests/TestAssertHandler.cpp \
	UnitTest++/tests/TestCheckMacros.cpp \
	UnitTest++/tests/TestChecks.cpp \
	UnitTest++/tests/TestCompositeTestReporter.cpp \
	UnitTest++/tests/TestCurrentTest.cpp \
	UnitTest++/tests/TestDeferredTestReporter.cpp \
	UnitTest++/tests/TestExceptions.cpp \
	UnitTest++/tests/TestLongMacros.cpp \
	UnitTest++/tests/TestMemoryOutStream.cpp \
	UnitTest++/tests/TestRequireMacrosWithExceptionsOff.cpp \
	UnitTest++/tests/TestRequireMacrosWithExceptionsOn.cpp \
	UnitTest++/tests/TestTest.cpp \
	UnitTest++/tests/TestTestList.cpp \
	UnitTest++/tests/TestTestMacros.cpp \
	UnitTest++/tests/TestTestResults.cpp \
	UnitTest++/tests/TestTestRunner.cpp \
	UnitTest++/tests/TestTestSuite.cpp \
	UnitTest++/tests/TestTimeConstraint.cpp \
	UnitTest++/tests/TestTimeConstraintMacro.cpp \
	UnitTest++/tests/TestUnitTestPP.cpp \
	UnitTest++/tests/TestXmlTestReporter.cpp

objects = $(patsubst %.cpp, %.o, $(src))
objects64 = $(patsubst UnitTest++/%, lib64/%, $(objects))
objects32 = $(patsubst UnitTest++/%, lib32/%, $(objects))
test_objects = $(patsubst %.cpp, %.o, $(test_src))
dependencies64 = $(subst .o,.d,$(objects64))
dependencies32 = $(subst .o,.d,$(objects32))
test_dependencies = $(subst .o,.d,$(test_objects))

define make-depend
  $(CCDEP) $(CCDEPFLAGS) -M $1 | \
  $(SED) -e 's,\($(notdir $2)\) *:,$(dir $2)\1: ,' > $3.tmp
  $(SED) -e 's/#.*//' \
      -e 's/^[^:]*: *//' \
      -e 's/ *\\$$//' \
      -e '/^$$/ d' \
      -e 's/$$/ :/' $3.tmp >> $3.tmp
  $(MV) $3.tmp $3
endef


all: $(test)

lib32: $(lib32)

$(lib64): $(objects64) 
	@echo Creating $(lib64) library...
	ar cr $(lib64) $(objects64)
	$(RANLIB) $(lib64)

$(lib32): $(objects32) 
	@echo Creating $(lib32) library...
	ar cr $(lib32) $(objects32)
	$(RANLIB) $(lib32)

$(test): $(lib64) $(test_objects)
	@echo Linking $(test)...
	$(CXX) $(LDFLAGS) -o $(test) $(test_objects) $(lib64) $(LDLIBS)
	@echo Running unit tests...
	./$(test)

clean:
	-@$(RM) -f $(objects64) $(objects32) $(test_objects) $(lib64) $(lib32) 2> /dev/null
	-@$(RM) -f $(dependencies64) $(dependencies32) $(test_dependencies) $(test) 2> /dev/null

lib64/%.o : UnitTest++/%.cpp
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CCFLAGS) -c -o $@ $<

lib32/%.o : UnitTest++/%.cpp
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CCFLAGS) -m32 -c -o $@ $<

%.o : %.cpp
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	$(CXX) $(CCFLAGS) -c -o $@ $<


ifneq "$(MAKECMDGOALS)" "clean"
-include $(dependencies64)
-include $(dependencies32)
-include $(test_dependencies)
endif
