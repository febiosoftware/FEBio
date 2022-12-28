import os
from os import path

OSNAME = "WINDOWS"

NUMCORES = os.cpu_count()

BUILDCOMMAND = ["MSBuild", "/p:configuration=Release"]
CLEANCOMMAND = ["MSBuild", "/p:configuration=Release", "/t:clean"]

FEBIODIR = os.getcwd()
FEBIOBUILDDIR = path.join(FEBIODIR, "cmbuild")
FEBIOPATH = path.join(FEBIOBUILDDIR, "bin", "Release", "febio4.exe")
FEBIOUPLOADPATH = FEBIOPATH
# FEBIOREMOTEDIR = "Windows/stage/febio"

# FBSDIR = "C:\\Users\\steve\\source\\repos\\cmbuild\\FEBioStudio\\"
# FBSBUILDDIR = FBSDIR + "build\\"
# FBSPATH = FBSBUILDDIR + "bin\\Release\\FEBioStudio.exe"
# FBSREMOTEDIR = "Windows/stage/bin"
# 
# CHEMDIR = "Z:\\FEBioChem\\"
# CHEMBUILDDIR = CHEMDIR + "cmbuild\\"
# CHEMPATH = CHEMBUILDDIR + "Release\\FEBioChem.dll"
# CHEMREMOTEDIR = "Windows/stage/febio"
# 
# HEATDIR = "Z:\\FEBioHeat\\"
# HEATBUILDDIR = HEATDIR + "cmbuild\\"
# HEATPATH = HEATBUILDDIR + "Release\\FEBioHeat.dll"
# HEATREMOTEDIR = "Windows/stage/febio"
# 
TESTDIR = path.join(FEBIODIR, "TestSuite/")
VERIFYDIR = path.join(TESTDIR, "Verify3/")
LOGDIR = path.join(TESTDIR, "Logs")
WORKINGDIR = VERIFYDIR

AUTOMATIONDIR = path.join(FEBIODIR, "ci/test-suite/")
GOLDSTANDARDS = "windowsGoldStandards.py"
RELEASELOG = path.join(AUTOMATIONDIR, "release.log")


# REMOTESCRIPT = "/root/update2/FEBioStudioDev/makeDevReleaseWindows.sh"

# pluginPaths = {'heat': HEATPATH, 'chem': CHEMPATH}

localExemptTests = []
