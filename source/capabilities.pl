#!/usr/bin/perl
# this script checks whether the installed g++ supports certain features we need
$res = "";
$os = `uname -s`;
$version = `$ARGV[0] --version 2> /dev/null`;
# remove comments in parentheses, they may or may not contain spaces
$version =~ s/\(.*?\)//g;
@v1 = split( /\s+/, $version );
$compiler = $v1[0];
# this is needed for CC
$compiler =~ s/://;
# g++ may identify itself as e.g. g++-7
$compiler =~ s/^g\+\+.*/g++/;
# at this point $compiler will hold "g++", "clang", "icc", or "CC" for supported compilers
if( $compiler eq "g++" ) {
    @v2 = split( /\./, $v1[1] );
#     $v2[0] is major revision, $v2[1] is minor revision, $v2[2] is microversion
#     precompiling headers is supported from g++ 3.4.0 onwards
#     it's broken on FreeBSD 11.0-RELEASE, so disable it on all BSD variants to be safe
#     it is also broken when using the openMPI mpiCC wrapper, so test for that as well
    if( $os !~ /bsd/i && $ARGV[0] !~ /mpi/ && (( $v2[0] == 3 && $v2[1] >= 4 ) || $v2[0] >= 4) ) {
	$res .= "precompile ";
    }
#     on modern FreeBSD installations, LLVM is the standard compiler
#     if you use GNU g++, you need different runtime support libraries
    if( $os =~ /FreeBSD/ ) {
	$libdir = "/usr/local/lib/gcc";
	if( $v2[0] <= 4 ) {
	    $libdir .= "$v2[0]$v2[1]";
	}
	else {
	    $libdir .= "$v2[0]";
	}
	if( -d $libdir ) {
	    $res .= "rpath=$libdir ";
	}
    }
#     vectorization is supported from g++ 4.0.0 onwards
    if( $v2[0] >= 4 ) {
	$res .= "vectorize ";
    }
#     ignore AVX capabilities on Darwin since assembler doesn't support it
    if( $os !~ /Darwin/ ) {
#         -march=native is supported from g++ 4.2.0 onwards, however we only
#         use it from 4.6.0 onwards to avoid sub-optimal AVX compilations...
	if( ( $v2[0] == 4 && $v2[1] >= 6 ) || $v2[0] >= 5 ) {
	    $res .= "native ";
	}
    }
}
# remove trailing spaces
$res =~ s/ +$//;
if( $res ne "" ) {
    print "$res\n";
}
