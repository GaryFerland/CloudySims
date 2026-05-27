#!/bin/bash

tempnam='vhtest.tmp.R6sq9'
export LD_LIBRARY_PATH="$PWD/../lib64"

test_cks_file () {
	rm -f $tempnam
	eval "arr=($1)"
	"${arr[@]}" > $tempnam 2>&1
	local retval1=$?
	if [ $retval1 -ne 0 ]; then
	   echo "command ==$1== failed"
	   exit 1;
	fi
	diff -q $tempnam $2
	local retval2=$?
	if [ $retval2 -ne 0 ]; then
	   echo "verification of file $2 failed"
	   exit 1;
	fi
	rm -f $tempnam
}

test_cks_stdin () {
	local msg=`cat $2 | $1`
	local cks1=`echo $msg | awk '{print $1}'`
	local msg2=`grep $2 $3`
	local cks2=`echo $msg2 | awk '{print $1}'`
	if [ $cks1 != $cks2 ]; then
		echo "checksum mismatch for file $2: got: $cks1, expected: $cks2"
		exit 1;
	fi
}

check_cmd () {
	$1 > /dev/null
	local retval1=$?
	if [ $retval1 -ne 0 ]; then
	   echo "command ==$1== failed"
	   exit 1;
	fi
}

check_error_msg () {
	rm -f $tempnam
	eval "arr=($1)"
	"${arr[@]}" > $tempnam 2>&1
	local retval1=$?
	if [ $retval1 -ne 1 ]; then
	   echo "command ==$1== did not fail"
	   exit 1;
	fi
	grep "$2" $tempnam > /dev/null
	local retval2=$?
	if [ $retval2 -ne 0 ]; then
		echo "error msg does not match ==$2=="
		exit 1;
	fi
	rm -f $tempnam
}

check_no_error_msg () {
	rm -f $tempnam
	eval "arr=($1)"
	"${arr[@]}" > $tempnam 2>&1
	grep "$2" $tempnam > /dev/null
	local retval2=$?
	if [ $retval2 -eq 0 ]; then
		echo "error msg is matched ==$2=="
		exit 1;
	fi
	rm -f $tempnam
}

check_no_output () {
	rm -f $tempnam
	eval "arr=($1)"
	"${arr[@]}" > $tempnam 2>&1
	if [ -s $tempnam ]; then
		echo "output was found in $tempnam"
		exit 1;
	fi
	rm -f $tempnam
}

test_cks_file "../bin/vh256sum -l 32 -b test*" "output_32.txt"
test_cks_file "../bin/vh128sum -bl64 test*" "output_64.txt"
test_cks_file "../bin/vh128sum --bin test*" "output_128.txt"
test_cks_file "../bin/vh256sum --bina test*" "output_256.txt"
test_cks_file "../bin/vh128sum --binary --length 512 test*" "output_512.txt"
test_cks_file "../bin/vh512sum -l 1024 -b test*" "output_1024.txt"

test_cks_file "../bin/vh256sum -l 32 --tag test*" "BSD_output_32.txt"
test_cks_file "../bin/vh128sum --length 64 --tag test*" "BSD_output_64.txt"
test_cks_file "../bin/vh128sum --tag test*" "BSD_output_128.txt"
test_cks_file "../bin/vh128sum --tag -l256 test*" "BSD_output_256.txt"
test_cks_file "../bin/vh128sum --tag --length 512 -- test*" "BSD_output_512.txt"
test_cks_file "../bin/vh512sum -l 1024 --tag test*" "BSD_output_1024.txt"

test_cks_file "../bin/vh256sum -l 32 --tag --scalar test*" "BSD_output_32.txt"
test_cks_file "../bin/vh128sum --length 64 --scalar --tag test*" "BSD_output_64.txt"
test_cks_file "../bin/vh128sum --tag --scalar test*" "BSD_output_128.txt"
test_cks_file "../bin/vh128sum --tag --scalar -l256 test*" "BSD_output_256.txt"
test_cks_file "../bin/vh512sum --scalar --tag -- test*" "BSD_output_512.txt"
test_cks_file "../bin/vh512sum -l 1024 --scalar --tag test*" "BSD_output_1024.txt"

test_cks_file "../bin/vh128sum --zero --binary test*" "output_zero_128.txt"
test_cks_file "../bin/vh128sum --tag -z test*" "BSD_output_zero_128.txt"

test_cks_stdin "../bin/vh256sum -l 32 -b" "test0128" "output_32.txt"
test_cks_stdin "../bin/vh256sum -l 32 -b -" "test0256" "output_32.txt"
test_cks_stdin "../bin/vh256sum -l 32 -b --scalar" "test0256" "output_32.txt"
test_cks_stdin "../bin/vh256sum -l64 -b" "test0128" "output_64.txt"
test_cks_stdin "../bin/vh256sum -l64 -b -" "test0512" "output_64.txt"
test_cks_stdin "../bin/vh256sum -l64 -b --scalar" "test0512" "output_64.txt"
test_cks_stdin "../bin/vh128sum -b" "test0128" "output_128.txt"
test_cks_stdin "../bin/vh128sum -b -" "test1024" "output_128.txt"
test_cks_stdin "../bin/vh128sum -b --scalar" "test1024" "output_128.txt"
test_cks_stdin "../bin/vh256sum -b" "test0128" "output_256.txt"
test_cks_stdin "../bin/vh256sum -b -" "test1536" "output_256.txt"
test_cks_stdin "../bin/vh256sum -b --scalar" "test1536" "output_256.txt"
test_cks_stdin "../bin/vh512sum -b" "test0128" "output_512.txt"
test_cks_stdin "../bin/vh512sum -b -" "test2048" "output_512.txt"
test_cks_stdin "../bin/vh512sum -b --scalar" "test2048" "output_512.txt"
test_cks_stdin "../bin/vh512sum -l 1024 -b" "test0128" "output_1024.txt"
test_cks_stdin "../bin/vh512sum -l 1024 -b -" "test3072" "output_1024.txt"
test_cks_stdin "../bin/vh512sum -l 1024 -b --scalar" "test3072" "output_1024.txt"

check_cmd "../bin/vh256sum -l 32 -c output_32.txt"
check_cmd "../bin/vh128sum -l64 -c output_64.txt"
check_cmd "../bin/vh128sum -c output_128.txt"
check_cmd "../bin/vh128sum -l256 -c output_256.txt"
check_cmd "../bin/vh128sum --leng 512 --check output_512.txt"
check_cmd "../bin/vh512sum -l 1024 --check output_1024.txt"

check_cmd "../bin/vh256sum -l 32 -cQ BSD_output_32.txt"
check_cmd "../bin/vh128sum -cQ -l64 BSD_output_64.txt"
check_cmd "../bin/vh128sum -cq BSD_output_128.txt"
check_cmd "../bin/vh128sum -l256 -c --quiet BSD_output_256.txt"
check_cmd "../bin/vh128sum --length 512 -Q --check BSD_output_512.txt"
check_cmd "../bin/vh512sum -l 1024 --check --stat BSD_output_1024.txt"

check_cmd "../bin/vh256sum -h"
check_cmd "../bin/vh512sum --help"
check_cmd "../bin/vh256sum -v"
check_cmd "../bin/vh512sum --vers"
check_cmd "../bin/vh512sum --verb test0128"
check_cmd "../bin/vh512sum --verbose --scalar test0128"

# make sure the code continues after a missing file
check_error_msg "../bin/vh256sum -l64 tost0000 -- --arg '' test0000" ": tost0000: No such file or directory"
check_error_msg "../bin/vh256sum -l64 tost0000 -- --arg '' test0000" ": --arg: No such file or directory"
check_error_msg "../bin/vh256sum -l64 tost0000 -- --arg '' test0000" ": '': No such file or directory"
check_error_msg "../bin/vh256sum -l64 tost0000 -b -- --arg '' test0000" "73711a77d6031b6f \*test0000"

# test widths that are not a power of 2
check_error_msg "../bin/vh256sum -l96 t0 test0000" "fdcdf538d6031b6f8a613d7f  test0000"
check_error_msg "../bin/vh256sum -l160 t0 test0000" "a49af44fdf9952d380d2cb84aace19bf6c0f58c8  test0000"
check_error_msg "../bin/vh256sum -l192 t0 test0000" "5481034fd2a9074478684acf85dceecd94e0e7830f82775b  test0000"
check_error_msg "../bin/vh256sum -l224 t0 test0000" "cfd4bf9386cf00766add0e4a30fc6e7ca6cf10f922cdf8e6cc5578ac  test0000"

check_error_msg "../bin/vh256sum -l64 -G" "invalid option -- 'G'"
check_error_msg "../bin/vh256sum -l64 --ver" "unrecognized option '--ver'"
check_error_msg "../bin/vh256sum -l64 -l sum" "invalid length: 'sum'"
check_error_msg "../bin/vh256sum -l64 -l 156" "invalid length: '156'"
check_error_msg "../bin/vh256sum -l64 -l156" "length must be a multiple of 32 between 32 and 1024"
check_error_msg "../bin/vh256sum -l64 --len sum" "invalid length: 'sum'"
check_error_msg "../bin/vh256sum -l64 --length 156" "invalid length: '156'"
check_error_msg "../bin/vh256sum -l64 --length 156" "length must be a multiple of 32 between 32 and 1024"
check_error_msg "../bin/vh256sum -l64 tost0000" "tost0000: No such file or directory"

check_error_msg "../bin/vh256sum -l64 -i test0128" "the --ignore-missing option is meaningful only when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -Q test0128" "the --status option is meaningful only when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -q test0128" "the --quiet option is meaningful only when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -w test0128" "the --warn option is meaningful only when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -s test0128" "the --strict option is meaningful only when verifying checksums"
check_error_msg "../bin/vh256sum -l64 --tag -c BSD_output_64.txt" "the --tag option is meaningless when verifying checksums"
check_error_msg "../bin/vh256sum -l64 --tag -t test0128" ": --tag does not support --text mode"
check_error_msg "../bin/vh256sum -l64 -ct output_64.txt" "the --binary and --text options are meaningless when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -c --text output_64.txt" "the --binary and --text options are meaningless when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -cb output_64.txt" "the --binary and --text options are meaningless when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -c --binary output_64.txt" "the --binary and --text options are meaningless when verifying checksums"
check_error_msg "../bin/vh128sum -cz output_zero_128.txt" "the --zero option is not supported when verifying checksums"
check_error_msg "../bin/vh256sum -l64 -cQ --verbose output_64.txt" "the --verbose option conflicts with --status"

check_error_msg "../bin/vh256sum -l 32 --check --strict error0_32.txt" "WARNING: 1 line is improperly formatted"

check_error_msg "../bin/vh128sum --check error1_128.txt" "test0768: OK"
check_error_msg "../bin/vh128sum --check error1_128.txt" "test1024: FAILED"
check_error_msg "../bin/vh128sum --check error1_128.txt" "tost3072: No such file or directory"
check_error_msg "../bin/vh128sum --check error1_128.txt" "tost3072: FAILED open or read"
check_error_msg "../bin/vh128sum --check error1_128.txt" "WARNING: 1 listed file could not be read"
check_error_msg "../bin/vh128sum --check error1_128.txt" "WARNING: 1 computed checksum did NOT match"
check_error_msg "../bin/vh128sum --check error1_128.txt" "WARNING: 1 line is improperly formatted"

check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "test0256: OK"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "test0128: FAILED"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "test2048: FAILED"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "tost0000: No such file or directory"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "tost0000: FAILED open or read"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "tost1536: No such file or directory"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "tost1536: FAILED open or read"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "WARNING: 2 listed files could not be read"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "WARNING: 2 computed checksums did NOT match"
check_error_msg "../bin/vh128sum -cl256 error2_256.txt" "WARNING: 2 lines are improperly formatted"

check_error_msg "../bin/vh128sum -cl64 error3_64.txt" "error3_64.txt: no properly formatted VH64 checksum lines found"

check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "test0256: OK"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "tost1536: No such file or directory"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "tost1536: FAILED open or read"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "test3072: FAILED"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "WARNING: 1 listed file could not be read"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "WARNING: 1 computed checksum did NOT match"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "WARNING: 1 line is improperly formatted"

check_error_msg "../bin/vh128sum --check --warn error1_128.txt" "error1_128.txt: 11: improperly formatted VH128 checksum line"
check_error_msg "../bin/vh128sum -wcl256 error2_256.txt" "error2_256.txt: 12: improperly formatted VH256 checksum line"
check_error_msg "../bin/vh128sum -wcl256 error2_256.txt" "error2_256.txt: 13: improperly formatted VH256 checksum line"
check_error_msg "../bin/vh128sum -wcl512 BSD_error1_512.txt" "BSD_error1_512.txt: 7: improperly formatted VH512 checksum line"

check_no_error_msg "../bin/vh128sum --check --ignore-missing error1_128.txt" "tost3072: No such file or directory"
check_no_error_msg "../bin/vh128sum -ci error1_128.txt" "tost3072: FAILED open or read"

check_no_output "../bin/vh128sum -Qcl256 error2_256.txt"
check_no_output "../bin/vh128sum -cl64 --status error3_64.txt"
check_no_output "../bin/vh128sum -cq output_128.txt"

echo "===================================="
echo "$0: all tests succeeded"
echo "===================================="
