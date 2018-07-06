set -v

if [ $# == 0 ]; then
	echo "Usage: Mkdir.bash platform"
	exit
fi

mkdir $1
cd $1
mkdir FEBio2
mkdir FEBioLib
mkdir FEBioMech
mkdir FEBioMix
mkdir FEBioOpt
mkdir FEBioPlot
mkdir FEBioTest
mkdir FEBioFluid
mkdir FEBioXML
mkdir FECore
mkdir NumCore


