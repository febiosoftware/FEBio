set -v

if [ $# == 0 ]; then
	echo "Usage: Mkdir.bash platform"
	exit
fi

cd FEBio2
mkdir $1
cd ..
cd FEBioHeat
mkdir $1
cd ..
cd FEBioLib
mkdir $1
cd ..
cd FEBioMech
mkdir $1
cd ..
cd FEBioMix
mkdir $1
cd ..
cd FEBioOpt
mkdir $1
cd ..
cd FEBioPlot
mkdir $1
cd ..
cd FEBioXML
mkdir $1
cd ..
cd FECore
mkdir $1
cd ..
cd NumCore
mkdir $1
cd ..
