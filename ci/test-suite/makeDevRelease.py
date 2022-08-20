#!/usr/bin/env python3

import sys, datetime
from tools import *

def tryBuild(name, baseDir, buildDir, path, branch=None):
    upOut = update(name, baseDir, branch)
    if not upOut[0]:
        print("\n\nSTDOUT:\n\n" + str(upOut[1].stdout, "Utf-8") + "\n\nSTDERR:\n\n" + str(upOut[1].stderr, "Utf-8"))
        sys.exit()
    
    buildOut = build(name, baseDir, buildDir, path)
    if not buildOut[0]:
        print("\n\nSTDOUT:\n\n" + str(buildOut[1].stdout, "Utf-8") + "\n\nSTDERR:\n\n" + str(buildOut[1].stderr, "Utf-8"))
        sys.exit()
        
    return buildOut[1]

dev = True
branch = "develop"
remoteDir = REMOTE_DEV_DIR
doUpload = True

if '-r' in sys.argv:
    dev = False
    branch = "master"
    remoteDir = REMOTE_RELEASE_DIR
    
if '-n' in sys.argv:
    doUpload = False



newFEBio = tryBuild("FEBio", FEBIODIR, FEBIOBUILDDIR, FEBIOPATH, branch)
newFBS = tryBuild("FEBioStudio", FBSDIR, FBSBUILDDIR, FBSPATH, branch)
newChem = tryBuild("FEBioChem", CHEMDIR, CHEMBUILDDIR, CHEMPATH)
newHeat = tryBuild("FEBioHeat", HEATDIR, HEATBUILDDIR, HEATPATH)

if OSNAME == "MACOS":
    macPostBuild()

if "-all" in sys.argv:
    newFEBio = True
    newFBS = True
    newChem = True
    newHeat = True
    
if "-febio" in sys.argv:
    newFEBio = True
    
if "-fbs" in sys.argv:
    newFBS = True
    
if "-chem" in sys.argv:
    newChem = True
    
if "-heat" in sys.argv:
    newHeat = True


if doUpload:
    if newFEBio or not dev:
        if not upload("FEBio", FEBIOUPLOADPATH, remoteDir + FEBIOREMOTEDIR):
            sys.exit()

    if newFBS or not dev:
        if not upload("FEBioStudio", FBSPATH, remoteDir + FBSREMOTEDIR):
            sys.exit()

    if newChem or not dev:
        if not upload("FEBioChem", CHEMPATH, remoteDir + CHEMREMOTEDIR):
            sys.exit()
            
    if newHeat or not dev:
        if not upload("FEBioHeat", HEATPATH, remoteDir + HEATREMOTEDIR):
            sys.exit()
            
    if newFEBio or newFBS or newChem or newHeat:
        if dev:
            if makeRelease():
                print("Release successfully made")
                
                logLines = []
                try:
                    with open(RELEASELOG, "r") as f:
                        logLines = f.readlines()
                except:
                    pass
                    
                newLine = datetime.datetime.now().strftime("%b %d %H:%M:%S - ")
                if newFEBio:
                    newLine += "FEBio, "
                if newFBS:
                    newLine += "FEBio Studio, "
                if newChem:
                    newLine += "FEBioChem, "
                if newHeat:
                    newLine += "FEBioHeat "
                    
                newLine = newLine.strip()
                if newLine[-1] == ',':
                    newLine = newLine [:-1]
                
                newLine += "\n"
                
                logLines.insert(0, newLine)
                
                with open(RELEASELOG, "w") as f:
                    for line in logLines[:100]:
                        f.write(line)
                
        else:
            print("Files built, uploaded, and staged. You must manually make a release on the server.")
    else:
        print("Everything is up to date. No need to upload")
