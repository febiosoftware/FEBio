import os
from os import path

OSNAME = "LINUX"

NUMCORES = os.cpu_count()

BUILDCOMMAND = ["make", "-j" + str(NUMCORES)]
CLEANCOMMAND = ["make", "clean"]

FEBIODIR = os.getcwd()
FEBIOBUILDDIR = path.join(FEBIODIR, "cmbuild")
FEBIOPATH = path.join(FEBIOBUILDDIR, "bin", "febio4")
FEBIOUPLOADPATH = FEBIOPATH

# FEBIOREMOTEDIR = "Linux/stage/bin"

# FBSDIR = "/home/sci/mherron/Projects/FEBioStudio/"
# FBSBUILDDIR = FBSDIR + "build/"
# FBSPATH = FBSBUILDDIR + "bin/FEBioStudio"
# FBSREMOTEDIR = "Linux/stage/bin"

# CHEMDIR = "/home/sci/mherron/Projects/Plugins/FEBioChem/"
# CHEMBUILDDIR = CHEMDIR + "cbuild/"
# CHEMPATH = CHEMBUILDDIR + "lib/libFEBioChem.so"
# CHEMREMOTEDIR = "Linux/stage/lib"
#
# HEATDIR = "/home/sci/mherron/Projects/Plugins/FEBioHeat/"
# HEATBUILDDIR = HEATDIR + "cbuild/"
# HEATPATH = HEATBUILDDIR + "lib/libFEBioHeat.so"
# HEATREMOTEDIR = "Linux/stage/lib"

TESTDIR = path.join(FEBIODIR, "TestSuite/")
VERIFYDIR = path.join(TESTDIR, "Verify4/")
LOGDIR = path.join(TESTDIR, "Logs")
WORKINGDIR = VERIFYDIR

AUTOMATIONDIR = path.join(FEBIODIR, "ci/test-suite/")
GOLDSTANDARDS = "linuxGoldStandards.py"
RELEASELOG = path.join(AUTOMATIONDIR, "release.log")

# REMOTESCRIPT = "/root/update2/FEBioStudioDev/makeDevReleaseLinux.sh"

# pluginPaths = {'heat': HEATPATH, 'chem': CHEMPATH}

localExemptTests = []
