/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

/* Configuration file specifically for Visual Studio -- every other system
   can generate this automatically so it's sure to be correct... */

#undef HAVE_POWI
#define HAVE_POW_DOUBLE_INT 1
#define HAVE_POW_DOUBLE_LONG 1
#define HAVE_POW_FLOAT_INT 1
#define HAVE_POW_FLOAT_LONG 1
#undef HAVE_POW_DOUBLE_FLOAT
#undef HAVE_POW_FLOAT_DOUBLE
#undef HAVE_FUNC
#define HAVE_ERF 1
#define HAVE_ERFC 1
#define HAVE_SINCOS 1
#define HAVE_LONG_LONG 1
#define HAVE_STRNLEN 1
#undef HAVE_CONSTEXPR
#undef HAVE_REALPATH
#undef HAVE_AVX_INTRIN
#undef HAVE_FMA_INTRIN
#undef HAVE_AVX2_INTRIN
#undef HAVE_AVX512F_INTRIN
