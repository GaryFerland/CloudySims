/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */
#include "cddefines.h"
#include "save.h"
#include "dynamics.h"
#include "grid.h"
#include "service.h"
#include "rfield.h"
#if defined(__unix) || defined(__APPLE__)
#include <cstddef>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#else
#define pid_t int
#define fork() TotalInsanityAsStub<pid_t>()
#define wait(X) TotalInsanityAsStub<pid_t>()
#endif
#ifdef __linux__
#include <sys/vfs.h>
// these come from /usr/include/linux/magic.h, but from linux kernel 3.7
// onwards, the file seems to have moved to /usr/include/uapi/linux/magic.h...
// to avoid the mess we copy the constants here
static const long NFS_SUPER_MAGIC = 0x6969L;
static const long SMB_SUPER_MAGIC = 0x517BL;
static const long AUTOFS_SUPER_MAGIC = 0x0187L;
// this one is not included in my header file (sshfs is also based on fuse)
static const long FUSE_SUPER_MAGIC = 0x65735546L;
#endif

STATIC void GridGatherOutput(const string&,long,const vector<long>&);
STATIC void GridGatherOutputSequential(const string&,long);
STATIC void GridGatherOutputParallel(const string&,long,const vector<long>&);
STATIC bool lgIsRemote(const string&);
STATIC void check_grid_file(const string&,int,int);

#ifndef MPI_ENABLED

int MPI_MODE_RDONLY = 0;
int MPI_MODE_WRONLY = 0;
int MPI_MODE_CREATE = 0;
int MPI_MODE_APPEND = 0;

int MPI_SUCCESS = 0;
int MPI_ERR_INTERN = -1;
MPI_File MPI_FILE_NULL = 0;

int total_insanity(MPI_File, int, MPI_Status*)
{
	return TotalInsanityAsStub<int>();
}

#endif

int mpi_mode_r = MPI_MODE_RDONLY;
int mpi_mode_w = MPI_MODE_WRONLY | MPI_MODE_CREATE;
int mpi_mode_a = MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND;

// NB NB this routine cannot throw any exceptions as it is executed outside
// the try{} block -- this includes mechanisms like ASSERT and cdEXIT!
void load_balance::init( unsigned int nJobs, unsigned int nCPU )
{
	if( nJobs == 0 )
		return;

	bool lgMPI = cpu.i().lgMPI();

	p_jobs.resize( nJobs );
	p_rank = cpu.i().nRANK();
	p_ptr = p_rank;
	p_ncpu = min(nCPU,nJobs);
	// the master rank will now set up a random sequence for the jobs
	// this way we hope to get statistical load balancing of the ranks
	if( p_rank == 0 )
	{
		for( unsigned int i=0; i < nJobs; ++i )
			p_jobs[i] = i;

		if( p_ncpu > 1 )
		{
			// This may or may not seed the random number generator used by random_shuffle.
			// The C++11 solution using <random> only works with g++ 6 or higher...
			srand( unsigned( time(NULL) ) );
			random_shuffle( p_jobs.begin(), p_jobs.end() );
		}
	}
	// now broadcast the random sequence to the other ranks...
	if( lgMPI )
		MPI_Bcast( &p_jobs[0], nJobs, MPI_type(p_jobs[0]), 0, MPI_COMM_WORLD );
	else
		for( unsigned int i=1; i < p_ncpu; ++i )
		{
			fflush(NULL);
			pid_t pid = fork();
			if( pid < 0 )
			{
				fprintf( ioQQQ, "creating the child process failed\n" );
				// this _should_ be exit() not cdEXIT()!
				exit(EXIT_FAILURE);
			}
			else if( pid == 0 )
			{
				/* this is child process */
				p_rank = i;
				p_ptr = p_rank;
				cpu.i().set_nRANK( p_rank );
				return;
			}
		}
}

void load_balance::finalize( exit_type exit_status )
{
	// wait for all jobs to finish
	if( cpu.i().lgMPI() )
		MPI_Barrier( MPI_COMM_WORLD );
	else
	{
		if( p_rank == 0 )
			for( unsigned int i=1; i < p_ncpu; ++i )
				(void)wait(NULL);
		else
			// this _should_ be exit() not cdEXIT()!
			exit(exit_status);
	}
}

/** process_output: concatenate output files produced in grid run */
void process_output()
{
	DEBUG_ENTRY( "process_output()" );

	// NOTE: when this routine is called all file handles have already been closed

	try
	{
		string main_output = save.chRedirectPrefix + ".out";

		// balance work over the ranks
		// rank n will process file numbers i with bound[n] <= i < bound[n+1]
		// in non-MPI runs, this will result in:
		//    bound[0] = 0;   bound[1] = grid.totNumModels;
		int nCPU = cpu.i().lgMPI() ? cpu.i().nCPU() : 1;
		int nRANK = cpu.i().lgMPI() ? cpu.i().nRANK() : 0;
		long stride = grid.totNumModels/nCPU;
		vector<long> bound( nCPU+1, stride );
		long remainder = grid.totNumModels - stride*nCPU;
		for( long i=1; i <= remainder; ++i )
			++bound[i];
		bound[0] = 0;
		for( long i=1; i <= nCPU; ++i )
			bound[i] += bound[i-1];

		ASSERT( bound[nCPU] == grid.totNumModels );

		// first process main output files
		if( !grid.lgKeepMainOutputSeparate )
			GridGatherOutput( main_output, grid.totNumModels, bound );

		// remove input files for individual grid points
		for( long j=0; j < grid.totNumModels; ++j )
		{
			if( j >= bound[nRANK] && j < bound[nRANK+1] )
			{
				string in_name = GridPointPrefix(j) + save.chRedirectPrefix + ".in";
				remove( in_name.c_str() );
			}
		}

		ASSERT( long(save.chFileName.size()) == save.nsave );

		for( long ipPun=0; ipPun < save.nsave; ++ipPun )
		{
			string fnam = save.chFilenamePrefix;
			fnam += save.chFileName[ipPun];
			// first do a minimal check on the validity of the save files
			for( int j=0; j < grid.totNumModels; ++j )
				if( j >= bound[nRANK] && j < bound[nRANK+1] )
					check_grid_file( fnam, j, ipPun );
			// and concatenate the output if necessary
			if( save.lgXSPEC[ipPun] )
			{
				if( cpu.i().lgMaster() )
				{
					ASSERT( save.FITStype[ipPun] >= 0 &&
						save.FITStype[ipPun] < NUM_OUTPUT_TYPES );

					// combine the grid.Spectra data from all ranks.
					// this is done by reading the results from file.
					for( int j=0; j < grid.totNumModels; ++j )
					{
						string gridnam = GridPointPrefix(j) + fnam;
						rd_block(&grid.Spectra[save.FITStype[ipPun]][j][0],
							 size_t(rfield.nflux)*sizeof(realnum),
							 gridnam.c_str());
						remove( gridnam.c_str() );
					}

					FILE *dest = open_data( fnam.c_str(), "ab" );
					// dest points to an empty file, so generate the complete FITS file now
					saveFITSfile( dest, save.FITStype[ipPun], save.punarg[ipPun][0],
						      save.punarg[ipPun][1], save.punarg[ipPun][2] );
					fseek( dest, 0, SEEK_END );
					ASSERT( ftell(dest)%2880 == 0 );
					fclose( dest );
				}
			}
			else if( save.lgSaveToSeparateFiles[ipPun] )
			{
				if( cpu.i().lgMaster() )
				{
					// open in binary mode in case we are writing a FITS file
					FILE *dest = open_data( fnam.c_str(), "ab" );
					// keep the save files for each grid point separate
					// the main save file contains the save header
					// salvage it by prepending it to the first save file
					string gridnam = GridPointPrefix(0) + fnam;
					append_file( dest, gridnam.c_str() );
					fclose( dest );
					// this will overwrite the old file gridnam
					rename( fnam.c_str(), gridnam.c_str() );
				}
			}
			else
			{
				GridGatherOutput( fnam, grid.totNumModels, bound );
			}
		}
	}
	catch( ... )
	{
		fprintf( ioQQQ, "PROBLEM - an internal error occurred while post-processing the grid output\n" );
	}
}

STATIC void GridGatherOutput(const string& basenam,
			     long nfiles,
			     const vector<long>& bound)
{
	if( cpu.i().lgMPI() && !lgIsRemote(basenam) )
		GridGatherOutputParallel( basenam, nfiles, bound );
	else {
		if( cpu.i().lgMaster() )
			GridGatherOutputSequential( basenam, nfiles );
	}
}

STATIC void GridGatherOutputSequential(const string& basenam,
				       long nfiles)
{
	// open in binary mode in case we are writing a FITS file
	FILE* output_handle = open_data( basenam.c_str(), "ab", AS_LOCAL_ONLY_TRY );
	if( output_handle == NULL )
		return;

	for( long j=0; j < nfiles; ++j )
	{
		string gridnam = GridPointPrefix(j) + basenam;
		append_file( output_handle, gridnam.c_str() );
		remove( gridnam.c_str() );
	}
	fclose( output_handle );
}

STATIC void GridGatherOutputParallel(const string& basenam,
				     long nfiles,
				     const vector<long>& bound)
{
	// determine total amount of data each rank has to copy
	// by summing the individual file sizes
	MPI_Offset mySize = 0;
	for( long j=0; j < nfiles; ++j )
	{
		if( j >= bound[cpu.i().nRANK()] && j < bound[cpu.i().nRANK()+1] )
		{
			string gridnam = GridPointPrefix(j) + basenam;
			FILE* fh = open_data( gridnam.c_str(), "r", AS_LOCAL_ONLY_TRY );
			if( fh != NULL )
			{
				fseek( fh, 0, SEEK_END );
				mySize += static_cast<MPI_Offset>( ftell(fh) );
				fclose(fh);
			}
		}
	}

	// broadcast the computed amounts to all ranks so that each
	// rank can compute the offset where it needs to start writing
	vector<MPI_Offset> offset(cpu.i().nCPU());
	for( int i=0; i < cpu.i().nCPU(); ++i )
	{
		MPI_Offset myCopy = mySize;
		// directly using &offset[i] below instead of the detour via
		// &myCopy leads to segfaults for reasons that I cannot fathom...
		MPI_Bcast( &myCopy, 1, MPI_type(myCopy), i, MPI_COMM_WORLD );
		offset[i] = myCopy;
	}

	MPI_File output_handle = open_data( basenam.c_str(), mpi_mode_a, AS_LOCAL_ONLY_TRY );
	if( output_handle == MPI_FILE_NULL )
		return;

	// compute offset where each rank needs to start writing
	MPI_Offset totalSize = 0;
	(void)MPI_File_get_size( output_handle, &totalSize );
	for( int j=0; j < cpu.i().nCPU(); ++j )
	{
		MPI_Offset tmp = offset[j];
		offset[j] = totalSize;
		totalSize += tmp;
	}

	// now gather the output and remove the individual files
	(void)MPI_File_set_view( output_handle, offset[cpu.i().nRANK()],
				 MPI_CHAR, MPI_CHAR, const_cast<char*>("native"),
				 MPI_INFO_NULL );
	for( long j=0; j < nfiles; ++j )
	{
		if( j >= bound[cpu.i().nRANK()] && j < bound[cpu.i().nRANK()+1] )
		{
			string gridnam = GridPointPrefix(j) + basenam;
			append_file( output_handle, gridnam.c_str() );
			remove( gridnam.c_str() );
		}
	}
	MPI_File_close( &output_handle );
}

// determine if a file resides on a remote share
// this is needed to determine if MPI-IO can operate safely on the file
#ifdef __linux__
STATIC bool lgIsRemote(const string& fnam)
{
	struct statfs buf;
	int res = statfs( fnam.c_str(), &buf );
	if( res != 0 )
		return true;
	// parallel file systems do not count as remote since MPI-IO is supported on those
	if( buf.f_type == NFS_SUPER_MAGIC ||
	    buf.f_type == SMB_SUPER_MAGIC ||
	    buf.f_type == AUTOFS_SUPER_MAGIC ||
	    buf.f_type == FUSE_SUPER_MAGIC )
		return true;
	else
		return false;
}
#else
STATIC bool lgIsRemote(const string&)
{
	// we do not know how to determine this, so we assume the worst
	return true;
}
#endif

/** check_grid_file: check whether the save file is present and is terminated by a GRID_DELIMIT string */
STATIC void check_grid_file( const string& fnam, int j, int ipPun )
{
	DEBUG_ENTRY( "check_grid_file()" );

	// these are binary files, don't touch them...
	if( save.lgFITS[ipPun] )
		return;

	bool lgForceNoDelimiter = false;
	// in these cases there should not be a GRID_DELIMIT string...
	if( !save.lgHashEndIter[ipPun] || !save.lg_separate_iterations[ipPun] ||
	    dynamics.lgTimeDependentStatic || strcmp( save.chHashString, "TIME_DEP" ) == 0 ||
	    strcmp( save.chHashString, "\n" ) == 0 )
		lgForceNoDelimiter = true;

	bool lgAppendDelimiter = true;
	bool lgAppendNewline = false;
	string gridnam = GridPointPrefix(j) + fnam;
	fstream str;
	open_data( str, gridnam.c_str(), mode_r, AS_LOCAL_ONLY_TRY );
	if( str.is_open() )
	{
		str.seekg( 0, ios_base::end );
		if( str.good() && str.tellg() > 0 )
		{
			// check if the file ends in a newline
			str.seekg( -1, ios_base::cur );
			char chr;
			str.get( chr );
			lgAppendNewline = ( chr != '\n' );
			// check if the GRID_DELIMIT string is present
			string line;
			str.seekg( 0, ios_base::beg );
			while( getline( str, line ) )
			{
				if( line.find( "GRID_DELIMIT" ) != string::npos )
					lgAppendDelimiter = false;
			}
		}
		str.close();
	}
	if( lgForceNoDelimiter )
		lgAppendDelimiter = false;
	if( lgAppendNewline || lgAppendDelimiter )
	{
		open_data( str, gridnam.c_str(), mode_a, AS_LOCAL_ONLY_TRY );
		if( str.is_open() )
		{
			if( lgAppendNewline )
				str << endl;
			if( lgAppendDelimiter )
			{
				str << save.chHashString << " GRID_DELIMIT -- grid";
				str << setfill( '0' ) << setw(9) << j << endl;
			}
			str.close();
		}
	}
}

/** append_file: append output produced on file <source> to open file descriptor <dest> */
void append_file( FILE *dest, const char *source )
{
	DEBUG_ENTRY( "append_file()" );

	FILE *src = open_data( source, "rb", AS_LOCAL_ONLY_TRY );
	if( src == NULL )
		return;

	// limited testing shows that using a 4 KiB buffer should
	// give performance that is at least very close to optimal
	// tests were done by concatenating 10 copies of a 62.7 MiB file
	const size_t BUF_SIZE = 4096;
	char buf[BUF_SIZE];

	while( ! feof(src) )
	{
		size_t nb = fread( buf, sizeof(char), BUF_SIZE, src );
		fwrite( buf, sizeof(char), nb, dest );
	}
	fclose(src);
	return;
}

/** append_file: append output produced on file <source> to open file handle <dest> using MPI I/O */
void append_file( MPI_File dest, const char *source )
{
	DEBUG_ENTRY( "append_file()" );

	FILE *src = open_data( source, "rb", AS_LOCAL_ONLY_TRY );
	if( src == NULL )
		return;

	// use larger buffer for parallel file systems
	const size_t BUF_SIZE = 32768;
	char buf[BUF_SIZE];

	while( ! feof(src) )
	{
		size_t nb = fread( buf, sizeof(char), BUF_SIZE, src );
		MPI_Status status;
		(void)MPI_File_write( dest, buf, nb, MPI_CHAR, &status );
	}
	fclose(src);
	return;
}
