#!/bin/bash
# Script to extract dependencies for makefiles

 
for dir in FEBio FECore; do
	cd $dir
	echo "" > dependencies.mk
	printf "OBJ = " > objects.mk
	for cppfile in $(ls *.cpp); do
		printf "${cppfile%%.*}.o : $cppfile" >> dependencies.mk
		printf "\\" >> objects.mk
		printf "\n${cppfile%%.*}.o " >> objects.mk
		for depend in $(grep \#include $cppfile | cut -f2 -d' '); do
			if [ ${depend:0:1} = \" ]; then
				depend=${depend:1}
				depend=${depend%%.*}
				if [ ${depend:0:7} = FECore/ ]; then
					depend="../"$depend
				fi
				if [ ${cppfile%%.*} = $depend ]; then
					printf " ${depend}.h" >> dependencies.mk
					hfile="${cppfile%%.*}.h"
					for dependh in $(grep \#include $hfile | cut -f2 -d' '); do
						if [ ${dependh:0:1} = \" ]; then
							dependh=${dependh:1}
							dependh=${dependh%%.*}
							if [ ${dependh:0:7} = FECore/ ]; then
								dependh="../"$dependh
							fi
							if [ -e "${dependh}.cpp" ]; then
								printf " ${dependh}.o" >> dependencies.mk
							else
								printf " ${dependh}.h" >> dependencies.mk
							fi
						fi
					done
				else
					if [ -e "${depend}.cpp" ]; then
						printf " ${depend}.o" >> dependencies.mk
					else
						printf " ${depend}.h" >> dependencies.mk
					fi
				fi
			fi
		done
		printf "\n" >> dependencies.mk
	done
#	cat dependencies.mk
	cd ..
done
cd FEBio
sed -i 's/ mat3d.o//; s/ windows.h//; s/ Console.h//' dependencies.mk
cd ../FECore
sed -i 's/ slu_ddefs.h pdsp_defs.h//; s/ windows.h//' dependencies.mk
exit
