# renames .F90 to .F95 in the current directory only

import os

SRCEXT=".F95"
DSTEXT=".F90"

print "Renaming files in",os.path.realpath(os.curdir)
files = os.listdir(os.curdir)
for file in files:
    root,ext = os.path.splitext( file )
    if( ext == SRCEXT ):
        os.rename( file, (root+DSTEXT) )
