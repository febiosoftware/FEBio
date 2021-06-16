/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include <FECore/fecore_api.h>
#include "febiomix_api.h"
#include <FECore/vec3d.h>

class FESolute;
class FEOsmoticCoefficient;
class FEMaterialPoint;

//------------------------------------------------------------------------
// This class should be used by all materials that support solutes
// TODO: This is a work in progress. The goal is to reduce the dynamic_casts to materials
//       that support solutes, and instead provide a single consistent interface to features
//       that need access to solute data (e.g. plot variables).
class FEBIOMIX_API FESoluteInterface
{
public:
	FESoluteInterface(){}
	virtual ~FESoluteInterface(){}

// derived classes need to implement the following functions
public:
	// return the number of solutes in the material
	virtual int Solutes() = 0;

	// return a solute material
	virtual FESolute* GetSolute(int i) = 0;

	// return the actual solution concentration at this material point
	virtual double GetActualSoluteConcentration(FEMaterialPoint& mp, int soluteIndex) { return 0.0; }

	// return the partition coefficient at this material point
	virtual double GetPartitionCoefficient(FEMaterialPoint& mp, int soluteIndex) { return 0.0; }

	// return the solute flux at this material point
	virtual vec3d GetSoluteFlux(FEMaterialPoint& mp, int soluteIndex) { return vec3d(0,0,0); }

	// get the osmotic coefficient
	virtual FEOsmoticCoefficient* GetOsmoticCoefficient() { return nullptr; }

	// get the osmolarity
	virtual double GetOsmolarity(const FEMaterialPoint& mp) { return 0.0; }

	// get the electric potential
	virtual double GetElectricPotential(const FEMaterialPoint& mp) { return 0.0; }

	// get the current density
	virtual vec3d GetCurrentDensity(const FEMaterialPoint& mp) { return vec3d(0,0,0); }

	// get the fixed charge density
	virtual double GetFixedChargeDensity(const FEMaterialPoint& mp) { return 0.0; }

	// get the referential fixed charge density
	virtual double GetReferentialFixedChargeDensity(const FEMaterialPoint& mp) { return 0.0; }

// additional member functions
public:
	// return the local index of a global solute ID (or -1, if the solute is not in this material)
	int FindLocalSoluteID(int soluteID);
};
