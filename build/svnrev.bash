#!/bin/bash

# Executing this bash script before compiling FEBio will add the svn revision number to the version
# number in the FEBio splash screen.  Note that this will only work if FEBio is under svn version control.

temp=$(svnversion -n)
svnrev=$( echo $temp | cut -f1 -d'M')

sed s/[$]WCREV[$]/$svnrev/ ../VS2010/svnrev_tpl.txt > ../FEBioLib/svnrev.h
