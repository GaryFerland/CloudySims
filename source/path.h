/* This file is part of Cloudy and is copyright (C)1978-2017 by Gary J. Ferland and
 * others.  For conditions of distribution and use see copyright notice in license.txt */

#ifndef PATH_H_
#define PATH_H_


/**\file path.h
 *
 * Cloudy needs to read a large number of data files when it starts.
 * The location of these data on your system is specified by 
 *
 * #define CLOUDY_DATA_PATH "XXXX" 
 *
 * below, so the code is set up so that data will be automatically
 * located.  "XXXX" gives the locations of the data files on your computer.
 *
 * This can also be set by using an argument to the compiler.  That is
 * discussed on the web site at 
 * http://wiki.nublado.org/wiki/EditPath
 *
 * The data path can also be set by the CLOUDY_DATA_PATH environment
 * variable, but it is generally easier to provide the permanent
 * location at compile time by editing CLOUDY_DATA_PATH below.
 */

/* The value below will be superseded if
 * specified by a compiler argument */
#ifndef CLOUDY_DATA_PATH 

/* 
 * Specify a search path to data directories here.
 *
 * This is for Unix variants.  You can specify a search path for one or several
 * data directories using the standard colon-separated list of Unix paths. You
 * can enter as many paths as you like. One example could be: */

//#define CLOUDY_DATA_PATH "/usr/local/cloudy/data:/home/user/cloudy/data"


/* The following (commented-out) example shows the format for Windows. You can
 * specify a search path for one or several data directories using a list separated
 * by semi-colon signs ";". The path must be enclosed between double quotes.
 * NB - note that the backslash "\" always needs to be typed twice, as shown below: */

//#define CLOUDY_DATA_PATH "c:\\projects\\cloudy\\trunk\\data;c:\\users\\gary\\data"

/* The // makes the above line of code a comment
 * if you want to use the Windows version remove the // from the line above and 
 * remove the unix version 12 lines above this line (or comment it out)
 * NB - make sure that the hash-sign "#" is in the first column */

#endif

/* That should be all you need to change! */

#endif /* PATH_H_ */
