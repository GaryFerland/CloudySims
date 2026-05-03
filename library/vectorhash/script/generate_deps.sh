#!/bin/sh

make_deps () {
	out=`echo $2 | sed s/:.*//`
	if [ "$1" != "src/vectorhash.cc" ]; then
		counter="${counter}="
		if [ "$counter" == "===" ]; then
			lib64="${lib64} \\\\\\n"
			lib32="${lib32} \\\\\\n"
			counter=""
		fi
		lib64="${lib64} lib64/$out"
		lib32="${lib32} lib32/$out"
	fi
	printf "lib64/$2\n"
	printf "\t\$(CXX) \$(CXXFLAGS) -c $1 -o \$@\n"
	printf "\n"
	printf "lib32/$2\n"
	printf "\t\$(CXX) \$(CXXFLAGS) -c -m32 $1 -o \$@\n"
	printf "\n"
}

make_deps_sub () {
	flag=""
	if [ $cpu = "intel" ] ; then
		case "$2" in
			*sse2*) flag="-msse2" ;;
			*avx2*) flag="-mavx2" ;;
			*avx512*) flag="-mavx512f" ;;
		esac
	fi
	out=`echo $3 | sed s/:.*//`
	counter="${counter}="
	if [ "$counter" == "===" ]; then
		lib64="${lib64} \\\\\\n"
		lib32="${lib32} \\\\\\n"
		counter=""
	fi
	lib64="${lib64} lib64/$1/$out"
	lib32="${lib32} lib32/$1/$out"
	printf "lib64/$1/$3\n"
	printf "\t\$(CXX) \$(CXXFLAGS) -c -D$1 $flag $2 -o \$@\n"
	printf "\n"
	printf "lib32/$1/$3\n"
	printf "\t\$(CXX) \$(CXXFLAGS) -c -m32 -D$1 $flag $2 -o \$@\n"
	printf "\n"
}

hardware () {
	case "$1" in
		x86*|amd64|i?86*|i86pc)
			cpu="intel"
			;;
		*)
			cpu="other"
			;;
	esac
}

installdirs () {
	libdir64="lib"
	libdir32="lib"
	if [ -d /usr/lib64 ] ; then
		libdir64="lib64"
	fi
	if [ -d /usr/lib32 ] ; then
		libdir32="lib32"
	fi
}

cxx=$1
cxxflags=$2

ver=1

OS=`uname -s`
hardware `uname -m`
installdirs

if [ "${OS}" = "Darwin" ] ; then
	ext="dylib"
	sflag="-x"
	nm64="libvhsum.${ver}.${ext}"
	nm32="libvhsum.${ver}.${ext}"
	libflags64="-dynamiclib -install_name @rpath/$nm64"
	libflags32="-dynamiclib -install_name @rpath/$nm32"
else
	ext="so"
	sflag=""
	nm64="libvhsum.${ext}.${ver}"
	nm32="libvhsum.${ext}.${ver}"
	libflags64="-shared"
	libflags32="-shared"
fi

lib64="lib64/${nm64}:"
lib32="lib32/${nm32}:"

counter=""

printf "LIB64DIR = $libdir64\n"
printf "\n"
printf "LIB32DIR = $libdir32\n"
printf "\n"
printf "LIB64NAME = lib64/$nm64\n"
printf "\n"
printf "LIB32NAME = lib32/$nm32\n"
printf "\n"
if [ "${OS}" = "Darwin" ] ; then
	printf "LDFLAGS = -lvhsum -Llib64 -Wl,-rpath,@executable_path/../\${LIB64DIR}\n"
else
	printf "LDFLAGS = -l:$nm64 -Llib64\n"
fi
printf "\n"
printf "lib64: \$(LIB64NAME)\n"
printf "\n"
printf "lib32: \$(LIB32NAME)\n"
printf "\n"
for file in src/*.cc
do
	multibuild=`grep -l 'EXT(' $file`
	deps=`$cxx $cxxflags -MM $file`
	if [ -n "$multibuild" ]; then
		make_deps_sub "VH32" "$file" "$deps"
		make_deps_sub "VH64" "$file" "$deps"
		make_deps_sub "VH128" "$file" "$deps"
		make_deps_sub "VH256" "$file" "$deps"
		make_deps_sub "VH512" "$file" "$deps"
		make_deps_sub "VH1024" "$file" "$deps"
	else
		make_deps "$file" "$deps"
	fi
done
for file in cpuid/*.cpp
do
	deps=`$cxx $cxxflags -MM $file`
	make_deps "$file" "$deps"
done
printf "$lib64\n"
printf "\t\$(CXX) \$(CXXFLAGS) $libflags64 -o lib64/$nm64 \$^\n"
printf "\tln -sf $nm64 lib64/libvhsum.${ext}\n"
printf "\n"
printf "$lib32\n"
printf "\t\$(CXX) \$(CXXFLAGS) -m32 $libflags32 -o lib32/$nm32 \$^\n"
printf "\tln -sf $nm32 lib32/libvhsum.${ext}\n"
printf "\n"
printf "install-lib:\n"
printf "\tmkdir -p \$(INSTALLDIR)/\$(LIB64DIR)\n"
printf "\tcp -af lib64/libvhsum.* \$(INSTALLDIR)/\$(LIB64DIR)\n"
printf "\tcd \$(INSTALLDIR)/\$(LIB64DIR); \\\\\n"
printf "\tstrip $sflag $nm64\n"
printf "\tmkdir -p \$(INSTALLDIR)/\$(LIB32DIR)\n"
printf "\tcp -af lib32/libvhsum.* \$(INSTALLDIR)/\$(LIB32DIR) 2> /dev/null || :\n"
printf "\tcd \$(INSTALLDIR)/\$(LIB32DIR); \\\\\n"
printf "\tstrip $sflag $nm32 2> /dev/null || :\n"
