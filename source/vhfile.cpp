/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#if ( defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__) ) && !defined(__CYGWIN__)
#include <unistd.h>
#else
#define _POSIX_MAPPED_FILES 0
#endif
#if _POSIX_MAPPED_FILES > 0
#include <sys/mman.h>
#endif

#include "vectorhash.h"
#include "cddefines.h"
#include "service.h"
#include "vhfile.h"

// width of the hash (in bits)
static const size_t vh_hash_width = 128;
// width of the hash in uint32_t elements
static const size_t vh_nstate = vh_hash_width/32;

// take checksum of a string
void VHstring(const string& s, void* out)
{
	DEBUG_ENTRY( "VHstring()" );

	VectorHash(s.data(), s.length(), 0xfd4c799d, out, vh_hash_width);
}

string VHstring(const string& s)
{
	DEBUG_ENTRY( "VHstring()" );

	uint32 cksum[vh_nstate];
	VHstring(s, cksum);

	ostringstream hash;
	for( uint32 i=0; i < vh_nstate; ++i )
		hash << hex << setfill('0') << setw(8) << cksum[i];

	return hash.str();
}

// take checksum of a file by supplying the full path
string VHfile(const string& fpath)
{
	DEBUG_ENTRY( "VHfile()" );

	auto fsize = FileSize(fpath);
	if( fsize == FS_UNKNOWN )
		return string();

	FILE* io = sys_fopen( fpath.c_str(), "r" );
	if( io == nullptr )
		return string();
#if _POSIX_MAPPED_FILES > 0
	int fd = fileno(io);
	char* map = (char*)mmap( nullptr, fsize, PROT_READ, MAP_SHARED, fd, 0 );
	if( fsize > 0 && map == MAP_FAILED )
		return fclose(io), string();
	uint32_t state[vh_nstate];
	VectorHash( map, fsize, 0xfd4c799d, state, vh_hash_width );
	munmap(map, fsize);
#else
	void* map = nullptr;
	if( fsize > 0 )
	{
		if( posix_memalign( &map, vh_hwreg_width/8, fsize ) != 0 )
			return fclose(io), string();
		if( fread( map, fsize, 1, io ) != 1 )
			return posix_memalign_free(map), fclose(io), string();
	}
	uint32_t state[vh_nstate];
	VectorHash( map, fsize, 0xfd4c799d, state, vh_hash_width );
	if( map != nullptr )
		posix_memalign_free(map);
#endif

	ostringstream hash;
	for( uint32_t i=0; i < vh_nstate; ++i )
		hash << hex << setfill('0') << setw(8) << state[i];

	fclose(io);
	return hash.str();
}

// this routine returns the checksum of a datafile. It filters out the eol characters,
// which makes it incompatible with the routine VHfile(), but also makes it OS agnostic...
// comment sections of lines starting with the hash symbol are also skipped
// this version is much slower than VHfile(), so should only be used on small files
// this routine uses the regular search path to find the file
string VHdatafile(const string& fname, access_scheme scheme)
{
	DEBUG_ENTRY( "VHdatafile()" );

	fstream ioFile;
	open_data(ioFile, fname, mode_r, scheme);

	string line, content;
	while( getline(ioFile, line) )
	{
		auto p = line.find('#');
		if( p != string::npos )
			line.erase(p);
		content += line;
	}

	return VHstring(content);
}
