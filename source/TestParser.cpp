/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "lines.h"
#include "parser.h"
#include "prt.h"

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
	TEST(TestTWavl)
	{
		LineSave.sig_figs = 6;

		// first test printing wavelengths in vacuum
		prt.lgPrintLineAirWavelengths = false;

		t_wavl t;
		CHECK( fp_equal(t.wavlVac(), -1_r) );
		t = 5000_vac;
		CHECK( fp_equal(t.wavlVac(), 5000_r) );
		CHECK( t.sprt_wl() == "5000.00A" );
		t = 6564.522546600_vac;
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6564.52A" );
		t = 6562.709693357_air;
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6564.52A" );
		t = t_vac(6564.522546600_r);
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6564.52A" );
		t = t_air(6562.709693357_r);
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6564.52A" );

		// now test printing wavelengths in air
		prt.lgPrintLineAirWavelengths = true;

		t = 6564.522546600_vac;
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6562.71A" );
		t = 6562.709693357_air;
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6562.71A" );
		t = t_vac(6564.522546600_r);
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6562.71A" );
		t = t_air(6562.709693357_r);
		CHECK( fp_equal(t.wavlVac(), 6564.522546600_r) );
		CHECK( t.sprt_wl() == "6562.71A" );

		// unary minus
		t = -t_vac(6564.522546600_r);
		CHECK( fp_equal(t.wavlVac(), -6564.522546600_r) );
	}
	TEST(TestReadWave)
	{
		// first test wavelengths in vacuum
		prt.lgPrintLineAirWavelengths = false;

		Parser p;
		p.setline("2316.23A");
		t_wavl t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.23_r) );
		CHECK( !p.lgEOL() );
		p.setline("12.1623m");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 12.1623e4_r) );
		p.setline("12.1623c");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 12.1623e8_r) );
		p.setline("2316.23 air");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.9419_r) );
		p.setline("2316.23 vacuum");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.23_r) );

		// now test wavelengths in air
		prt.lgPrintLineAirWavelengths = true;

		p.setline("2316.23A");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.9419_r) );
		p.setline("2316.23 air");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.9419_r) );
		p.setline("2316.23 vacuum");
		t = p.getWave();
		CHECK( fp_equal(t.wavlVac(), 2316.23_r) );

		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		p.setline("no number");
		t = p.getWaveOpt();
		CHECK( p.lgEOL() );
		p.setline("no number");
		CHECK_THROW( (void)p.getWave(), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}
	TEST(TestReadRange)
	{
		// first test printing wavelengths in vacuum
		prt.lgPrintLineAirWavelengths = false;

		Parser p;
		t_wavl t1, t2;
		p.setline("range 1200 2.3c");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 1200_r) );
		CHECK( fp_equal(t2.wavlVac(), 2.3e8_r) );
		p.setline("range 4100 air 5600 vacuum");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 4101.1571_r) );
		CHECK( fp_equal(t2.wavlVac(), 5600_r) );
		p.setline("range 4100 air 5600");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 4101.1571_r) );
		CHECK( fp_equal(t2.wavlVac(), 5600_r) );

		// now test wavelengths in air
		prt.lgPrintLineAirWavelengths = true;

		p.setline("range 1200 5400a");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 1200_r) );
		CHECK( fp_equal(t2.wavlVac(), 5401.5013_r) );
		p.setline("range 4100 air 5600 vacuum");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 4101.1571_r) );
		CHECK( fp_equal(t2.wavlVac(), 5600_r) );
		p.setline("range 4100 5600 vacuum");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( fp_equal(t1.wavlVac(), 4101.1571_r) );
		CHECK( fp_equal(t2.wavlVac(), 5600_r) );

		// test failure modes
		p.setline("1200 2.3c");
		CHECK( !p.GetRange("RANG", t1, t2) );
		p.setline("range");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( p.lgEOL() );
		p.setline("range 1200 vacuum");
		CHECK( p.GetRange("RANG", t1, t2) );
		CHECK( p.lgEOL() );
	}
	TEST(TestReadLineID)
	{
		// LineID always stores wavelength in vacuum internally, but what it reads from the line
		// can be either air or vacuum wavelength, depending on the value of prt.lgPrintLineAirWavelengths
		// explicit keywords AIR or VACUUM can also be used to force the interpretation of the wavelength
		prt.lgPrintLineAirWavelengths = false;

		Parser p;
		LineID line;
		p.setline("H  1  1216");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && line.wavlVac() == 1216_r && line.indLo() < 0 && line.indHi() < 0 );
		p.setline("H  1 1216A");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && line.wavlVac() == 1216_r && line.ELo() < 0_r );
		p.setline("\"Fe 2b\"  12.00m");
		line = p.getLineID();
		CHECK( line.chLabel() == "Fe 2b" && line.wavlVac() == 12.00e4_r );
		p.setline("\"Fe 2b  \"  12.00c # comment");
		line = p.getLineID();
		CHECK( line.chLabel() == "Fe 2b" && line.wavlVac() == 12.00e8_r );
		p.setline("CO  12.00C # comment");
		line = p.getLineID();
		CHECK( line.chLabel() == "CO" && line.wavlVac() == 12.00e8_r );
		p.setline("Al 2 1670. index=1,5");
		line = p.getLineID();
		CHECK( line.chLabel() == "Al 2" && line.wavlVac()== 1670_r && line.indLo() == 1 && line.indHi() == 5 && line.ELo() < 0_r );
		p.setline("Al 2 1670. Elow=1");
		line = p.getLineID();
		CHECK( line.chLabel() == "Al 2" && line.wavlVac() == 1670_r && line.indLo() < 0 && line.indHi() < 0 && line.ELo() == 1_r );
		p.setline("monitor line \"H  1\" 1216 index=1,5 -100.");
		line = p.getLineID(false);
		CHECK( line.chLabel() == "H  1" && line.wavlVac() == 1216_r && line.indLo() == 1 && line.indHi() == 5 );
		double x = p.FFmtRead();
		CHECK( x == -100. );
		p.setline("monitor line \"He 2\" 1217 -110.");
		line = p.getLineID(false);
		CHECK( line.chLabel() == "He 2" && line.wavlVac() == 1217_r && line.indLo() == -1 );
		x = p.FFmtRead();
		CHECK( x == -110. );

		p.setline("H  1 6562.71A air");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && fp_equal(line.wavlVac(), 6564.523_r) );
		CHECK( line.str() == "\"H  1\" 6564.52A" );
		p.setline("H  1 6562.71A air index=2, 5");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && fp_equal(line.wavlVac(), 6564.523_r) );
		CHECK( line.str() == "\"H  1\" 6564.52A" );
		CHECK( line.indLo() == 2 && line.indHi() == 5 );

		prt.lgPrintLineAirWavelengths = true;
		p.setline("H  1 6562.71A");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && fp_equal(line.wavlVac(), 6564.523_r) );
		p.setline("H  1 6564.52A vacuum");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && fp_equal(line.wavlVac(), 6564.52_r) );
		CHECK( line.str() == "\"H  1\" 6562.71A" );
		// test an input line with two line IDs
		p.setline("stop line \"c  2\" 157.636m air relative to \"o  3\" 5008.24 vacuum");
		line = p.getLineID(false);
		CHECK( line.chLabel() == "c  2" && fp_equal(line.wavlVac(), 157.67897e4_r) );
		CHECK( line.str() == "\"c  2\" 157.636m" );
		t_wavl t = line.twav();
		CHECK( fp_equal(line.wavlVac(), t.wavlVac()) );
		line = p.getLineID(false);
		CHECK( line.chLabel() == "o  3" && fp_equal(line.wavlVac(), 5008.24_r) );
		CHECK( line.str() == "\"o  3\" 5006.84A" );

		prt.lgPrintLineAirWavelengths = false;
		// test failure modes
		FILE *bak = ioQQQ;
		FILE *tmp = tmpfile();
		if( tmp != NULL )
			ioQQQ = tmp;
		// not really a failure, but a warning
		p.setline("H 1   1216a");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && line.wavlVac() == 1216._r );
		p.setline("\"H 1\"   1216M");
		line = p.getLineID();
		CHECK( line.chLabel() == "H  1" && line.wavlVac() == 1216.e4_r );
		// the rest are all real errors
		p.setline("monitor line \"H  1\" 1216");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("H 1");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("H  1");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("TOTL12.00m");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2b 12.00m");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("\"Fe 2b  12.00m");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  comment  12.00m");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  index");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  index=3");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  index=-1,2");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  index=3,2");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  elow");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("Fe 2  12.00m  elow=-1");
		CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		p.setline("normalize to Fe 2  12.00m");
		CHECK_THROW( (void)p.getLineID(false), cloudy_exit );
		//p.setline("Fe 2  12.00m  keyword");
		//CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		//p.setline("Fe 2  12.00 m");
		//CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		//p.setline("H  1 6564.52A index=2,3 vacuum");
		//CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		//p.setline("H  1 6564.52A Elow=82258.92 vacuum");
		//CHECK_THROW( (void)p.getLineID(), cloudy_exit );
		if( tmp != NULL )
			fclose(tmp);
		ioQQQ = bak;
	}
	TEST(TestReadInputTime)
	{
		Parser p;
		p.setline("age 1 second");
		CHECK_EQUAL(1.,parse_input_time(p));
		p.setline("age 1 minute");
		CHECK_EQUAL(60.,parse_input_time(p));
		p.setline("age 1 hour");
		CHECK_EQUAL(3600.,parse_input_time(p));
		p.setline("age 1 hour log");
		CHECK_EQUAL(36000.,parse_input_time(p));
		p.setline("age 1 days");
		CHECK_EQUAL(86400.,parse_input_time(p));
		p.setline("age 1 weeks");
		CHECK_EQUAL(604800.,parse_input_time(p));
		p.setline("age 2 fortnights");
		CHECK_EQUAL(2419200.,parse_input_time(p));
		p.setline("age 1 month");
		CHECK_EQUAL(2629800.,parse_input_time(p));
		p.setline("age 1 year");
		CHECK_EQUAL(31557600.,parse_input_time(p));
		p.setline("age 1 century");
		realnum age = 3155760000.;
		CHECK(fp_equal_tol(age,parse_input_time(p),1e-6_r*age));
		p.setline("age 1 millenia");
		age = 31557600000.;
		CHECK(fp_equal_tol(age,parse_input_time(p),1e-6_r*age));

		p.setline("stop time 1 Myr");
		age = 31557600000000.;
		CHECK(fp_equal_tol(age,parse_input_time(p),1e-6_r*age));
		p.setline("stop time 1 Gyr");
		age = 31557600000000000.;
		CHECK(fp_equal_tol(age,parse_input_time(p),1e-6_r*age));
	}
}
