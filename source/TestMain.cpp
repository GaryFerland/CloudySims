/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <fenv.h>
#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include "cddefines.h"

int main ()
{
	ioQQQ = stdout;
	// disable FP exceptions, they would lead to spurious crashes in the tests
	fesetenv(FE_DFL_ENV);
	return UnitTest::RunAllTests();
}
