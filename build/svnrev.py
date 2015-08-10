#! /usr/bin/env python
# -*- coding: utf-8 -*-
import subprocess, sys, os

# Use subwcrev on Windows to determine the svn revision and create svnrev.h

# Find the revision
os.chdir("..")
feb_dir = os.getcwd()
revision = subprocess.Popen(['svnversion', '-n'], stdout=subprocess.PIPE).communicate()[0]
if ':' in revision: revision = revision.split(":")[1]
if 'M' in revision: revision = revision.split("M")[0]

# Write svnrev.h
f_svnrev = open(feb_dir + "/FEBioLib/svnrev.h", "w")
f_svnrev.write("//  This file is created by svnrev.py which uses svnversion in Linux and OSX\n")
f_svnrev.write("//  to determine the svn revision number\n\n")
f_svnrev.write("#define SVNREVISION  " + revision + "\n")
f_svnrev.close()