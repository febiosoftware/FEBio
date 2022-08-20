OSNAME = "MACOS"

NUMCORES = 4

BUILDCOMMAND = ["make", "-j" + str(NUMCORES)]
CLEANCOMMAND = ["make", "clean"]

FBSDIR = "/Users/mherron/Projects/FEBioStudio/"
FBSBUILDDIR = FBSDIR + "build/"
FBSPATH = FBSBUILDDIR + "bin/FEBioStudio.app/Contents/MacOS/FEBioStudio"
FBSREMOTEDIR = "macOS/stage/FEBioStudio.app/Contents/MacOS"

FEBIODIR = "/Users/mherron/Projects/FEBio/"
FEBIOBUILDDIR = FEBIODIR + "cbuild/Release/"
FEBIOPATH = FEBIOBUILDDIR + "bin/febio3"
FEBIOUPLOADPATH = FBSBUILDDIR + "bin/FEBioStudio.app/Contents/MacOS/febio3"
FEBIOREMOTEDIR = "macOS/stage/FEBioStudio.app/Contents/MacOS"

CHEMDIR = "/Users/mherron/Projects/Plugins/FEBioChem/"
CHEMBUILDDIR = CHEMDIR + "cbuild/"
CHEMPATH = CHEMBUILDDIR + "lib/libFEBioChem.dylib"
CHEMREMOTEDIR = "macOS/stage/FEBioStudio.app/Contents/Frameworks"

HEATDIR = "/Users/mherron/Projects/Plugins/FEBioHeat/"
HEATBUILDDIR = HEATDIR + "cbuild/"
HEATPATH = HEATBUILDDIR + "lib/libFEBioHeat.dylib"
HEATREMOTEDIR = "macOS/stage/FEBioStudio.app/Contents/Frameworks"

TESTDIR = "/Users/mherron/Projects/TestSuite/"
VERIFYDIR = TESTDIR + "Verify3/"
LOGDIR = TESTDIR + "Logs/"
WORKINGDIR = "/Users/mherron/scratch/FEBioTests/"

AUTOMATIONDIR = "/Users/mherron/Projects/automation/"
GOLDSTANDARDS = "macOSGoldStandards.py"
RELEASELOG = AUTOMATIONDIR + "release.log"

REMOTESCRIPT = "/root/update2/FEBioStudioDev/makeDevReleaseMacOS.sh"

pluginPaths = {'heat': HEATPATH, 'chem': CHEMPATH}

localExemptTests = []
