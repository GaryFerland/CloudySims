/* This file is part of Cloudy and is copyright (C)1978-2025 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef VECTORHASH_H
#define VECTORHASH_H

#include <string>

// take checksum of aligned buffer
void VectorHash(const void* key, size_t len, uint32_t seed, void* out);
// take checksum of a file by supplying the full path
std::string VHfile(const std::string& fpath);

#endif
