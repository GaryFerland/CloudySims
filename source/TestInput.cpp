/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cdstd.h"
#include <UnitTest++.h>
#include "cddefines.h"
#include "input.h"

namespace {
	TEST(TestLgIsCommentSeq)
	{
		CHECK( !lgIsCommentSeq( "# comment", false, false ) );
		CHECK( lgIsCommentSeq( "## comment", false, false ) );
		CHECK( lgIsCommentSeq( "// comment", false, false ) );
		CHECK( lgIsCommentSeq( "% comment", false, false ) );
		CHECK( !lgIsCommentSeq( "* comment", false, false ) );
		CHECK( lgIsCommentSeq( "; comment", false, false ) );
		CHECK( !lgIsCommentSeq( "c omment", false, false ) );
		CHECK( !lgIsCommentSeq( "hden", false, false ) );

		CHECK( !lgIsCommentSeq( "# comment", true, false ) );
		CHECK( lgIsCommentSeq( "## comment", true, false ) );
		CHECK( lgIsCommentSeq( "// comment", true, false ) );
		CHECK( lgIsCommentSeq( "% comment", true, false ) );
		CHECK( lgIsCommentSeq( "* comment", true, false ) );
		CHECK( !lgIsCommentSeq( "; comment", true, false ) );
		CHECK( !lgIsCommentSeq( "c omment", true, false ) );
		CHECK( !lgIsCommentSeq( "hden", true, false ) );

		CHECK( lgIsCommentSeq( "# comment", false, true ) );
		CHECK( lgIsCommentSeq( "## comment", false, true ) );
		CHECK( lgIsCommentSeq( "// comment", false, true ) );
		CHECK( lgIsCommentSeq( "% comment", false, true ) );
		CHECK( !lgIsCommentSeq( "* comment", false, true ) );
		CHECK( lgIsCommentSeq( "; comment", false, true ) );
		CHECK( !lgIsCommentSeq( "c omment", false, true ) );
		CHECK( !lgIsCommentSeq( "hden", false, true ) );

		CHECK( lgIsCommentSeq( "# comment", true, true ) );
		CHECK( lgIsCommentSeq( "## comment", true, true ) );
		CHECK( lgIsCommentSeq( "// comment", true, true ) );
		CHECK( lgIsCommentSeq( "% comment", true, true ) );
		CHECK( lgIsCommentSeq( "* comment", true, true ) );
		CHECK( !lgIsCommentSeq( "; comment", true, true ) );
		CHECK( !lgIsCommentSeq( "c omment", true, true ) );
		CHECK( !lgIsCommentSeq( "hden", true, true ) );
	}

	TEST(TestLgInputComment)
	{
		CHECK( lgInputComment( "# comment" ) );
		CHECK( lgInputComment( "## comment" ) );
		CHECK( lgInputComment( "// comment" ) );
		CHECK( lgInputComment( "% comment" ) );
		CHECK( lgInputComment( "* comment" ) );
		CHECK( !lgInputComment( "; comment" ) );
		CHECK( !lgInputComment( "c omment" ) );
		CHECK( !lgInputComment( "hden" ) );
	}

	TEST(TestLgInputEOF)
	{
		CHECK( lgInputEOF( "" ) );
		CHECK( lgInputEOF( " text" ) );
		CHECK( lgInputEOF( "***" ) );
		CHECK( !lgInputEOF( "**" ) );
		CHECK( !lgInputEOF( "command" ) );
	}

	TEST(TestStripComment)
	{
		string line = "command # comment\n";
		StripComment( line, false );
		CHECK( line == "command # comment" );

		line = "command ## comment\n";
		StripComment( line, false );
		CHECK( line == "command " );

		line = "command // comment\n";
		StripComment( line, false );
		CHECK( line == "command " );

		line = "command % comment\n";
		StripComment( line, false );
		CHECK( line == "command " );

		line = "command * comment\n";
		StripComment( line, false );
		CHECK( line == "command * comment" );

		line = "command ; comment\n";
		StripComment( line, false );
		CHECK( line == "command " );

		line = "command c comment\n";
		StripComment( line, false );
		CHECK( line == "command c comment" );

		line = "command # comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command ## comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command // comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command % comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command * comment\n";
		StripComment( line, true );
		CHECK( line == "command * comment" );

		line = "command ; comment\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command c comment\n";
		StripComment( line, true );
		CHECK( line == "command c comment" );

		line = "command \"# ## // % * ; c text\"\n";
		StripComment( line, true );
		CHECK( line == "command \"# ## // % * ; c text\"" );

		line = "command #\"# ## // % * ; c text\"\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command \"# ## // % * ; c text\n";
		StripComment( line, true );
		CHECK( line == "command \"# ## // % * ; c text" );

		line = "command # \"# ## // % * ; c text\n";
		StripComment( line, true );
		CHECK( line == "command " );

		line = "command # \"# ## // % * ; c text\n";
		StripComment( line, false );
		CHECK( line == "command # \"# ## // % * ; c text" );

		input.lgUnderscoreFound = false;
		input.lgBracketFound = false;

		line = "command \"file_name[]\" _on_ # comment _ [ ]\n";
		StripComment( line, false );
		CHECK( line == "command \"file_name[]\"  on  # comment _ [ ]" );
		CHECK( input.lgUnderscoreFound && !input.lgBracketFound );

		line = "command param\r";
		StripComment( line, true );
		CHECK( line == "command param" );

		line = "# comment ## another comment\r";
		StripComment( line, false );
		CHECK( line == "# comment ## another comment" );

		line = "command # comment ## another comment\n";
		StripComment( line, false );
		CHECK( line == "command # comment ## another comment" );
	}

	TEST(TestGetString)
	{
		string line = "\"label \\ # ## // % * ; c text_23[]\"";
		string label;
		size_t p = GetString(line, 0, label);
		CHECK( label == "label \\ # ## // % * ; c text_23[]" );
		CHECK( p == line.length() );

		line = "\"text";
		p = GetString(line, 0, label);
		CHECK( label == "" );
		CHECK( p == string::npos );

		line = " \"text\"";
		CHECK_THROW( GetString(line, 0, label), bad_assert );
	}
}
