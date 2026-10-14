/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VHFILE_H
#define VHFILE_H

// take checksum of the contents of a string
void VHstring(const string& s, void* out);
string VHstring(const string& s);
// take checksum of a file by supplying the full path
string VHfile(const string& fpath);
// take checksum of a datafile by stripping eol chars and comment lines
string VHdatafile(const string& fnam, access_scheme scheme=AS_DEFAULT);

#endif
