/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "parser.h"

namespace {
	TEST(TestReadNumber)
	{
		Parser p;
		p.setline("1000");
		CHECK_EQUAL(1000,p.FFmtRead());
	}
	TEST(TestReadOffEnd)
	{
		Parser p;
		p.setline("1000");
		CHECK(!p.lgEOL());
		CHECK_EQUAL(1000,p.FFmtRead());
		CHECK(!p.lgEOL());
		CHECK_EQUAL(0,p.FFmtRead());
		CHECK(p.lgEOL());
	}
	TEST(TestReadNegativeNumber)
	{
		Parser p;
		p.setline("-1000");
		CHECK_EQUAL(-1000,p.FFmtRead());
	}
	TEST(TestReadPlusNumber)
	{
		Parser p;
		p.setline("+1000");
		CHECK_EQUAL(+1000,p.FFmtRead());
	}
	TEST(TestReadFraction)
	{
		Parser p;
		p.setline("0.125");
		CHECK_EQUAL(0.125,p.FFmtRead());
	}
	TEST(TestReadEmbeddedFraction)
	{
		Parser p;
		p.setline("Pi is 3.14159");
		CHECK(fp_equal_tol(3.14159,p.FFmtRead(),1e-3));
	}
	TEST(TestReadEmbeddedFractions)
	{
		Parser p;
		p.setline("Pi is 3.14159, e is 2.71828");
		CHECK(fp_equal_tol(3.14159,p.FFmtRead(),1e-3));
		CHECK(fp_equal_tol(2.71828,p.FFmtRead(),1e-3));
	}
//	TEST(TestReadCommaNumber)
//	{
//		Parser p;
//		p.setline("1,000");
//		CHECK_EQUAL(1000,p.FFmtRead());
//	}
	TEST(TestReadOKCommaNumber)
	{
		Parser p;
		p.setline("1000,");
		CHECK_EQUAL(1000,p.FFmtRead());
	}
	TEST(TestReadExponentialNumber)
	{
		Parser p;
		p.setline("1e3");
		CHECK_EQUAL(1000,p.FFmtRead());
	}
	TEST(TestReadExponentialFraction)
	{
		Parser p;
		p.setline("0.125e1");
		CHECK_EQUAL(1.25,p.FFmtRead());
	}
	TEST(TestReadPositiveExponentialFraction)
	{
		Parser p;
		p.setline("0.125e+1");
		CHECK_EQUAL(1.25,p.FFmtRead());
	}
	TEST(TestReadNegativeExponentialFraction)
	{
		Parser p;
		p.setline("1.25e-1");
		CHECK_EQUAL(0.125,p.FFmtRead());
	}
	TEST(TestReadExponentNumber)
	{
		Parser p;
		p.setline("10^3,");
		CHECK_EQUAL(1000,p.FFmtRead());
	}
	TEST(TestReadFractionalExponentNumber)
	{
		Parser p;
		p.setline("10.^3.5");
		CHECK(fp_equal_tol(3162.277,p.FFmtRead(),1e-2));
	}
	TEST(TestReadFractionalSquaredNumber)
	{
		Parser p;
		p.setline("2.5^2");
		CHECK_EQUAL(6.25,p.FFmtRead());
	}
	TEST(TestReadFractionalSquaredNegativeNumber)
	{
		// At present unary - binds tighter that the exponential
		Parser p;
		p.setline("-2.5e0^2e0");
		CHECK_EQUAL(6.25,p.FFmtRead());
	}
	TEST(TestReadChainedExponentNumber)
	{
		Parser p;
		p.setline("10^2^3");
		CHECK_EQUAL(1e8,p.FFmtRead());
	}
	TEST(TestReadProductNumber)
	{
		Parser p;
		p.setline("1.25*5.0");
		CHECK_EQUAL(6.25,p.FFmtRead());
	}
	TEST(TestReadProductPowExpr)
	{
		Parser p;
		p.setline("1.25*10^2");
		CHECK_EQUAL(125,p.FFmtRead());
	}
	TEST(TestReadPowProductExpr)
	{
		Parser p;
		p.setline("10^2*1.25");
		CHECK_EQUAL(125,p.FFmtRead());
	}
	TEST(TestReadProductProductExpr)
	{
		Parser p;
		p.setline("2*2*1.25");
		CHECK_EQUAL(5,p.FFmtRead());
	}
	TEST(TestReadDivExpr)
	{
		Parser p;
		p.setline("4/2");
		CHECK_EQUAL(2,p.FFmtRead());
	}
	TEST(TestReadDivDivExpr)
	{
		Parser p;
		p.setline("9/2/2");
		CHECK_EQUAL(2.25,p.FFmtRead());
	}
	TEST(TestReadDivMulExpr)
	{
		Parser p;
		p.setline("9/2*2");
		CHECK_EQUAL(9,p.FFmtRead());
	}
	TEST(TestReadMulDivExpr)
	{
		Parser p;
		p.setline("2*9/2");
		CHECK_EQUAL(9,p.FFmtRead());
	}
	TEST(TestReadExpDivExpExpr)
	{
		Parser p;
		p.setline("3^3/2^2");
		CHECK_EQUAL(6.75,p.FFmtRead());
	}
	TEST(TestReadProductPowProductExpr)
	{
		Parser p;
		p.setline("2*10^2*1.25");
		CHECK_EQUAL(250,p.FFmtRead());
	}
	TEST(TestReadMultiProductProductExpr)
	{
		Parser p;
		p.setline("3*10*10*10*10*10");
		CHECK_EQUAL(3e5,p.FFmtRead());
		p.setline("10*10*10*10*10*3");
		CHECK_EQUAL(3e5,p.FFmtRead());
	}
	TEST(TestReadVariable)
	{
		Parser p;
		p.setline("$a=5");
		p.doSetVar();
		p.setline("$col=6");
		p.doSetVar();
		p.setline("$a");
		CHECK_EQUAL(5,p.FFmtRead());
		p.setline("$a*$col");
		CHECK_EQUAL(5*6,p.FFmtRead());
		p.setline("$col*$a");
		CHECK_EQUAL(5*6,p.FFmtRead());
		p.setline("$a*7");
		CHECK_EQUAL(5*7,p.FFmtRead());
		p.setline("7*$a");
		CHECK_EQUAL(5*7,p.FFmtRead());
		p.setline("2^$a");
		CHECK_EQUAL(32,p.FFmtRead());
	}
}
