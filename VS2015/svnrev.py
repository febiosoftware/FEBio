#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess, sys, os

# Use subwcrev on Windows to determine the svn revision and create svnrev.h

# Find the revision
proj_dir = os.getcwd()
#print(proj_dir)
os.chdir("../..")
feb_dir = os.getcwd()
#print(feb_dir)
sub_out = subprocess.Popen(['subwcrev', feb_dir], stdout=subprocess.PIPE).communicate()[0]
#print(sub_out)
revision_line = sub_out.split("\n")[2]
if 'Mixed' in revision_line: revision = revision_line.split(":")[1]
elif 'Updated' in revision_line: revision = revision_line.split(" ")[3]
else: sys.exit("Error: unknown subwcrev output")
#print(revision)

# Write svnrev.h
f_svnrev = open(feb_dir + "/FEBioLib/svnrev.h", "w")
f_svnrev.write("//  This file is created by svnrev.py which uses subwcrev in Windows\n")
f_svnrev.write("//  to determine the svn revision number\n\n")
f_svnrev.write("#define SVNREVISION  " + revision + "\n")
f_svnrev.close()