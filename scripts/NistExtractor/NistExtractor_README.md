NistExtractor
=============
NistExtractor.py retrieves the energy levels and transition probabilities from NIST for a single species. 
The desired species is the first parameter.
The second and optional parameter allows the user to specify an early cutoff to limit the number of levels to be used by Cloudy



Example: Fe IX  with 30 levels

```
NistExtractor.py Fe_IX 30
```

NEWrapper.py allows for NistExtractor to be run on multiple species.
The required parameter is a file listing the desired species and optionally the level limit.

See "all_species.txt" and "cloudy_species.txt" as examples of how to format the input file.

```

# -----------------------------------------------------------------------------
Originally written by Matt Lykins in support of:
     Lykins et al. 2015, ApJ, 807, 118  [doi:10.1088/0004-637X/807/1/118]

 Updated and extended by Maryam Dehghanian (June 2025)in support of:
     JWST Archival Research Program AR-6019

 Description:
 This script constructs a STOUT-style atomic data directory from NIST ASD data.
 It now accepts a wide range of ion name formats, including:
     "O_III", "o_iii", and "o_3"
 All of the above are correctly interpreted as O$^{2+}$ (i.e., O+2).

Example on how to run the code:

```
Python3 NistExtractor.py O_III
```

 Notes:
 - The output directory is structured to be compatible with the Cloudy STOUT database.
 - Since NIST does not provide collision strengths, the generated .coll file is empty.
   If used in Cloudy, electron collisions will be estimated using the g-bar approximation.
# -----------------------------------------------------------------------------
