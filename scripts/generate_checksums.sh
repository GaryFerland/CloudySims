#!/bin/sh
echo -n "generating checksums.dat... "
rm -f checksums.dat
export LD_LIBRARY_PATH=`realpath ../library/vectorhash/lib64`
find . -type f -exec ../library/vectorhash/bin/vh128sum '{}' '+' > checksums.tmp6R4sQl
sed 's$  ./$  $' checksums.tmp6R4sQl | grep -v 'checksums.tmp6R4sQl' | LC_COLLATE=C sort -k 2 > checksums.dat
rm -f checksums.tmp6R4sQl
echo "done."
