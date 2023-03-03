#!/usr/bin/perl

# This perl script cleans out the test suite after a test run, it will
# only leave files that are part of the official distribution and the
# cloudy executable (with a name ending in .exe) if that was present!
#
# Peter van Hoof

while( defined( $input = glob("*") ) ) {
	@ll = split( /\./, "$input" );
	if( $#ll != 1 || ( $ll[1] ne "in" && $ll[1] ne "pl" && $ll[1] ne "pm" && $ll[1] ne "htm" && $ll[1] ne "jpg" && $ll[1] ne "dat" && $ll[1] ne "txt" && $ll[1] ne "exe" && $ll[1] ne "pdf" && $ll[1] ne "vsz" && $ll[1] ne "md" ) ) {
		print "deleting $input...\n";
		unlink "$input";
	}
	elsif( $ll[0] =~ /^grid[0-9]{9}/ ) {
		print "deleting $input...\n";
		unlink "$input";
	}
}
unlink "mpi_batch_jobs.txt";
unlink "checkend.txt";
unlink "close.txt";
unlink "crashed.txt";
unlink "debug.txt";
unlink "minor.txt";
unlink "serious.txt";
unlink "skip.txt";
unlink "wvlng.txt";
unlink "warnings.txt";
unlink "optimal.in";
