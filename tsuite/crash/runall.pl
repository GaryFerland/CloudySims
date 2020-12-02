#!/usr/bin/perl -w
# the above is generic form of #! for perl, with warnings
 
# of the cloudy executable on your system

# various VS execs
#$exe = "c:/Projects/Cloudy/trunk/debug/trunk.exe";
#$exe = "c:/projects/Cloudy/trunk/release/trunk.exe";

#$exe = "c:/Projects/Cloudy/trunk/vs08/debug/vs08.exe";
$exe = "c:/projects/cloudy/trunk/vs08/release/vs08.exe";

# icc 
#$exe = "c:/projects/cloudy/trunk/icl/release/cloudy_icl.exe";

$nMod = 0;
$nSkip = 0;

# this loops over all the *.in files in the current directory
# and runs the code to produce *.out files
while ( defined( $input = glob("*.in") ) )
{
	if ( $nMod>=$nSkip )
	{
                # $output =~ s/.out//gi;
                # $output =~ s/.out//gi;
		print( "$input going to " );
		$output = $input;
		$output =~ s/\.in//gi;
		if( -e "$output".".out" )
		{
			rename( "$output".".out" , "$output".".bak" );
		}
                $out = "$output".".out";
		print("$output\n");
		# actually execute the code
#		system "nice -n 5 $exe < $input  > $out";
		system "$exe < $input  > $out";
	}
	++$nMod;
}

printf( "\n %i models were computed, and %i were skipped.\n ",$nMod,$nSkip);

