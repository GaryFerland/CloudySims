/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"

namespace {

	struct LongIntFixture
	{
		auto_vec<long> p;
		LongIntFixture()
		{
			p = auto_vec<long>( new long[10] );
			for( int i=0; i < 10; ++i )
				p[i] = i;
		}
		~LongIntFixture() {}
		auto_vec<long> myfunc()
		{
			auto_vec<long> a( new long[10] );
			for( int i=0; i < 10; ++i )
				a[i] = 2*i;
			return a;
		}
	};

	struct MyClassFixture
	{
		struct Class
		{
			long n;
			Class() { n = 23; }
			~Class() {}
		};
		auto_vec<Class> p;
		MyClassFixture()
		{
			p = auto_vec<Class>( new Class[10] );
		}
		~MyClassFixture() {}
		auto_vec<Class> myfunc()
		{
			auto_vec<Class> a( new Class[10] );
			for( int i=0; i < 10; ++i )
				a[i].n = i+17;
			return a;
		}
		auto_vec<Class> myfunc2()
		{
			return auto_vec<Class>();
		}
	};

	TEST_FIXTURE(LongIntFixture,TestConstructBasic)
	{
		// test the basic constructor
		auto_vec<long> q;
		CHECK( q.get() == NULL );
		q = auto_vec<long>( new long[5] );

		// note that q is not used below and valgrind
		// should check whether the memory in q is freed

		// this tests assigning an auto_vec to another
		auto_vec<long> r( p );
		CHECK( p.get() == NULL );
		CHECK( r.data() != NULL );
		// alternative way of doing the same
		auto_vec<long> t;
		t = r;
		CHECK( r.get() == NULL );
		CHECK( t.data() != NULL );
		// finally, returning a vector as function result!
		t = myfunc();
		// this checks operator[]
		for( int i=0; i < 10; ++i )
			CHECK_EQUAL(2*i,t[i]);
		// construct straight from function result
		auto_vec<long> c( myfunc() );
		for( int i=0; i < 10; ++i )
			CHECK_EQUAL(2*i,c[i]);
	}

	// now repeat the tests with a class that has a constructor
	TEST_FIXTURE(MyClassFixture,TestConstructClass)
	{
		auto_vec<Class> q;
		CHECK( q.get() == NULL );
		q = auto_vec<Class>( new Class[5] );
		// check whether the constructor was executed
		for( int i=0; i < 5; ++i )
			CHECK_EQUAL(23,q[i].n);
		auto_vec<Class> r( p );
		CHECK( p.get() == NULL );
		CHECK( r.data() != NULL );
		auto_vec<Class> t;
		t = q;
		CHECK( q.get() == NULL );
		CHECK( t.data() != NULL );
		t = myfunc();
		for( int i=0; i < 10; ++i )
			CHECK_EQUAL(i+17,t[i].n);
		auto_vec<Class> a,b,c;
		a = b = c = t;
		CHECK( b.get() == NULL );
		CHECK( c.get() == NULL );
		CHECK( t.get() == NULL );
		for( int i=0; i < 10; ++i )
			CHECK_EQUAL(i+17,a[i].n);
		a.reset();
		CHECK( a.get() == NULL );
		a = myfunc2();
		CHECK( a.get() == NULL );
		a = auto_vec<Class>();
		CHECK( a.get() == NULL );
	}

}
