Two types of files are described here.  

The *.abn files are used to specify the chemical composition.
One of these files will be used if the command like the following appears:
abundances "filename.abn"
These abn files may also specify grains, although grains are not included by default.
default.abn specifies our default composition.  default-reference.abn is a copy of that file.

The *.dep files specify depletion factors that account for the elements that are used to build grains.
These may be parsed by the
metals deplete "filename.dep"
or 
metals deplete Jenkins 2009 Fstar=0.5 "filename.dep"
commands.

==================================================
default "solar" composition
default.abn gives the default composition used if no others are specified.
default-backup.abn is a backup of this file, use to reestablish default if you change it
This is a fairly old set of "solar" abundances that we maintain for backwards compatibility

To change the default abundance set replace the default.abn file with your 
own set of default abundances. Copy default-backup.abn to default.abn to reestablish 
our default abundances

Other .abn files are available - see the comments within each file for more details.
The file name gives a good indication of its contents

=======================================================
default depletion files

metals deplete
by default will read the ISM_CloudyClassic.dep files.  The format should be clear.
If the element is not specified then its depletion is set to unity.

metals deplete Jankins2009 Fstar
reads the ISM_Jenkins09_Tab4.dep file.  If this file is updated or replaced then
its current format must be maintained.

=======================================================


Format of the *.abn files

comments begin with #

The file parsing ends with a line containing a field of stars, as
****************
Everything after the field of stars is ignored, so remaining lines can be used
to document the abundances, in addition to comment lines starting with #.

The grains command can be included.  The abundance files recognize all grain command options.
This makes it possible to include grains in your mix of gas and dust.

Names of elements are followed by the abundance by number relative to hydrogen.

Hydrogen may be in specified, or may not be.  If the abundance of hydrogen is 
given and is not equal to unity the code will rescale the abundances by the entered value 
for hydrogen.  If hydrogen is not specified it is assumed to have an abundance
of 1.

All elements heavier than hydrogen are turned off before the abn files are read.
An element is turned on if it appears in the abn file.
If an element is not included in the *.abn file it will not be included in the calculation.

==========================================================

Original coding by Joshua Schlueter, NSF REU 2012 Summer
