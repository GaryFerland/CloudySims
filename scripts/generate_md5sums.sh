#!/bin/sh
echo -n "generating md5datafiles.dat... "
rm -f md5datafiles.dat
find . -type f -exec ../source/md5datafile.exe '{}' ';' > md5datafiles.tmp6R4sQl
sed 's$  ./$  $' md5datafiles.tmp6R4sQl | grep -v 'md5datafiles.tmp6R4sQl' > md5datafiles.dat
rm -f md5datafiles.tmp6R4sQl
echo "done."
