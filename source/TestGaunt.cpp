#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "physconst.h"
#include "atmdat_gaunt.h"

namespace {

	struct GauntFixture
	{
		double SanityCheckGaunt(long Z, double loggam2, double logu, double refval, double relerr)
		{
			double gam2 = exp10(loggam2);
			double u = exp10(logu);
			double ZZ = double(Z);
			double Te = TE1RYD*pow2(ZZ)/gam2;
			double anu = pow2(ZZ)*u/gam2;
			double val = t_gaunt::Inst().gauntff( Z, Te, anu );
			return fabs(val/refval-1.)/relerr;
		}
	};

	TEST_FIXTURE(GauntFixture,TestGaunt)
	{
		// our Gaunt factors are merged relativistic and non-relativistic data
		// for high gamma^2 we have non-relativistic data, these are compared to Sutherland (1998)
		// for low gamma^2 we have relativistic data, these are compared to Nozawa+ (1998)

		// these data were taken from
		// >>refer	HI	gff	Sutherland R.S., 1998, MNRAS, 300, 321
		CHECK( SanityCheckGaunt( 1, 4., -4., 2.7008, 0.0004 ) <= 1. );
		CHECK( SanityCheckGaunt( 1, 4., 4., 1.1040, 0.0004 ) <= 1. );
		CHECK( SanityCheckGaunt( 30, 4., 0., 1.0237, 0.0004 ) <= 1. );
		CHECK( SanityCheckGaunt( 12, 3., -3., 2.2126, 0.0004 ) <= 1. );
		CHECK( SanityCheckGaunt( 17, 1., -1., 1.5088, 0.0004 ) <= 1. );
		// these data were taken from
		// >>refer	gff	Nozawa S., Itoh N., Kohyama Y., 1998, ApJ, 507, 530
		CHECK( SanityCheckGaunt( 1, -4., -4., 6.120, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 2, -4., -4., 8.934, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 6, -4., -4., 28.92, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 7, -4., -4., 34.50, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 8, -4., -4., 40.18, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 1, -3., -1., 1.852, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 2, -3., -1., 1.921, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 6, -3., -1., 3.399, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 7, -3., -1., 4.050, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 8, -3., -1., 4.772, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 6, -2.5, 2., 107.7, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 7, -2.5, 2., 199.0, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 8, -2.5, 2., 330.8, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 1, 0., 3., 0.2001, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 2, 0., 3., 0.2041, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 6, 0., 3., 0.3340, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 7, 0., 3., 0.4358, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 8, 0., 3., 0.5909, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 10, -3., 2., 3925., 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 12, -3., 2., 6118., 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 14, -3., 2., 8653., 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 16, -3., 2., 11470., 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 26, -3., 2., 28710., 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 10, 1., 3., 0.5375, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 12, 1., 3., 0.5603, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 14, 1., 3., 0.5962, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 16, 1., 3., 0.6497, 0.005 ) <= 1. );
		CHECK( SanityCheckGaunt( 26, 1., 3., 1.370, 0.005 ) <= 1. );
	}
}
