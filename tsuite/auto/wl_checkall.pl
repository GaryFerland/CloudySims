#!/usr/bin/perl
#
#  wl_checkall.pl:
# 	Check all output files and write all unique wavelength misses per sim to
#  file.  Totals (including unique lines across all sims) are printed at stdout.
#

use strict;
use warnings;

my $outfile = "wvlng.txt";


#
#  Conversion factor for wavelength comparison.  This is required to match the
#  intended wavelength reported by Cloudy with the one that appears in the script
#  itself.
#  The factor is meant to be applied to the first of the input wavelengths.
#
sub convert_wl
{
	my ($wl_unit1, $wl_unit2) = @_;


	my $conv_fact = 1;
	if( $wl_unit2 ne $wl_unit1 )
	{
		if( $wl_unit2 eq "A" or $wl_unit2 eq "")
		{
			$conv_fact = 1e4	if( $wl_unit1 eq "m" );
			$conv_fact = 1e6	if( $wl_unit1 eq "c" );
		}
		elsif( $wl_unit2 eq "m" )
		{
			$conv_fact = 1e-4	if( $wl_unit1 eq "A" or $wl_unit1 eq "" );
			$conv_fact = 1e2	if( $wl_unit1 eq "c" );
		}
		elsif( $wl_unit2 eq "c" )
		{
			$conv_fact = 1e-6	if( $wl_unit1 eq "A" or $wl_unit1 eq "" );
			$conv_fact = 1e-2	if( $wl_unit1 eq "m" );
		}
	}
	#	print "$wl_unit2\t  VS  \t$unit_oldwl\t $conv_fact:\t $scrpt_wl\t VS ".	($oldwl * $conv_fact)	."\n";

	return	$conv_fact;
}



#
#  Process wavelength to report the numerical value, the unit, and the number of
#  significant figures found.
#
sub wl_proc
{
	my ($wl) = @_;

	my $wl_unit = "";
	$wl_unit = substr($wl, length($wl)-1, 1, "")	if ($wl =~ m/[a-zA-Z]$/);

	my $idot = index( $wl, "." );
	my $nsigfig = 0;
	if( $idot >= 0 )
	{
		$nsigfig = length(substr($wl, $idot+1));
	}

	return	($wl, $wl_unit, $nsigfig);
}



#
#  Given two wavelengths (including units), produce a string of the first
#  wavelength that best matches the format of the second wavelength.  This is
#  used to pattern match the wavelength reported by Cloudy against the one
#  that exists in the script (likely of lower precision).
#
sub wl_match
{
	my ($wl1, $wl2) = @_;

	my ($wl1str, $wl1_unit, $wl1_nsigfig) = &wl_proc($wl1);
	my ($wl2str, $wl2_unit, $wl2_nsigfig) = &wl_proc($wl2);
	my $conv_fact = &convert_wl($wl1_unit, $wl2_unit);

	return	sprintf("%.*f", $wl2_nsigfig, $wl1str * $conv_fact);
}



#
#  Report if two wavelengths match after taking into account differences in
#  units and precision.
#
sub do_wl_match
{
	my ($wl1, $wl2) = @_;

	my (undef, $wl2_unit) = &wl_proc($wl2);
	my $match_wl = &wl_match($wl1, $wl2);

	my $match = 0;
	$match = 1	if( $wl2 eq "$match_wl$wl2_unit" );
	#	print "$wl1 \t $wl2\t $match_wl\t $match\n";

	return	($match);
}



sub get_file_contents
{
	my ($script) = @_;

	open SCRIPT, "< $script" or die  "Could not open:\t $script\n";
	my @contents = <SCRIPT>;
	close SCRIPT	or print "Could not close:\t $script\n";

	return  @contents;
}



#
#  Produce errors for report.
#
sub prep_report_fails
{
	my ($output) = @_;

	my $script = $output;
	$script =~ s/\.out$/.in/;

	my @contents = &get_file_contents( $output );
	my @fix_list;

	# Uncomment for old reporting scheme of 3 lines of test per mismatch.
	#	my $FindLine_Fail_Mesg = "PROBLEM findline did not find line with label";
	#	my $FindLine_Suggestion = "The closest with correct label was";
	my $FindLine_Fail_Mesg = "WARNING: no exact matching lines found for";
	my $FindLine_Suggestion = "Taking best match as";

	for( my $i = 0; $i < scalar(@contents); $i++ )
	{
		my $line = $contents[$i];

		if( $line =~ $FindLine_Fail_Mesg )
		{
			chomp($line);
			$line =~ s/\.$//;

			my %fix;

			my @words = split(/\"/, $line);
			$fix{error} = $line;
			$fix{label} = $words[1];
			$fix{oldwl} = $words[2];
			$fix{oldwl} =~ s/\s+//g;
			#	@words      = split(/\s+/, $words[2]);
			#	$fix{oldwl} = $words[3];
			#	$fix{oldwl} =~ s/\.$//;
			#	print "label = ". $fix{label} ."\n";
			#	print "oldwl = ". $fix{oldwl} ."\n";

			# Uncomment for old reporting scheme of 3 lines of test per mismatch.
			#	++$i;
			$line = $contents[++$i];

			if( $line !~ $FindLine_Suggestion )
			{
				print	"Script: $script\t WARNING!  No closest match found! "
				  .	"Ignore line:\t" . $fix{label} ." @".$fix{oldwl} ."\n";
				next;
			}

			@words	    = split(/\"/, $line);
			@words	    = split(/\s+/, $words[1]);
			$fix{newwl} = $words[$#words];
			#	print "newwl = ". $fix{newwl} ."\n";

			# Needed for error checking...
			$fix{found} = 0;
			$fix{matches} = [];

			#	print "label= ". $fix{label} ."\t oldwl = ". $fix{oldwl} ."\t ". $fix{newwl} ."\n";
			
			my $is_listed = 0;
			for( my $i = 0; $i < scalar(@fix_list); $i++ )
			{
				if( $fix_list[$i]{label} eq $fix{label}  and
				    &do_wl_match($fix_list[$i]{oldwl}, $fix{oldwl}) and
				    &do_wl_match($fix_list[$i]{newwl}, $fix{newwl}) )
				{
					$is_listed++;
					last;
				}
			}

			if( not $is_listed )
			{
				$fix{error} .= "\t closest wavelength:\t $fix{newwl}\n";
				push( @fix_list, { %fix } );
			}
		}
	}

	return	\@fix_list;
}






open OUTFILE, "> $outfile"	or die "Error: Could not open: $outfile\n";


my @AllFixes;
my ($nlines, $nfiles, $nuniq) = (0, 0, 0);
foreach my $output ( glob "*.out" )
{
	my $fix_list = &prep_report_fails( $output );
	if( @$fix_list )
	{
		for( my $i = 0; $i < scalar(@$fix_list); $i++)
		{
			print OUTFILE "$output: ". $$fix_list[$i]{error};
		}

		$nlines += @$fix_list;
		$nfiles++;

		for (my $i = 0; $i < scalar(@$fix_list); $i++)
		{
			my $found = 0;
			for( my $ifix = 0; $ifix < scalar(@AllFixes); $ifix++ )
			{
				if( $AllFixes[$ifix]{label} eq $$fix_list[$i]{label} and
				    &do_wl_match($AllFixes[$ifix]{oldwl}, $$fix_list[$i]{oldwl} ) and
				    &do_wl_match($AllFixes[$ifix]{newwl}, $$fix_list[$i]{newwl} ) )
				{
					$found++;
					last;
				}
			}

			if( not $found )
			{
				my %fix;
				$fix{label} = $$fix_list[$i]{label};
				$fix{oldwl} = $$fix_list[$i]{oldwl};
				$fix{newwl} = $$fix_list[$i]{newwl};
				push( @AllFixes, { %fix } );
				$nuniq++;
			}
		}
	}
}

close OUTFILE	or warn "Warning: Could not close: $outfile\n";


print	"\nLooking for unmatched wavelengths:";
if( -s $outfile )
{
	print	"\nWARNING! In $nfiles sims: $nlines unidentified wavelengths ($nuniq unique across all sims).\n";
	print	  "WARNING! See file:  $outfile\n";
}
else
{
	print	" good, none found.\n";
}

print	"\n";
