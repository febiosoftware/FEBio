/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include <FECore/FEMaterial.h>
#include "FEFluid.h"
#include "FESolutesMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for FluidFSI materials.

class FEBIOFLUID_API FEFluidSolutesMaterial2 : public FEMaterial, public FESoluteInterface
{
public:
	FEFluidSolutesMaterial2(FEModel* pfem);

	// returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData() override;

public:
	FEFluid* GetFluidMaterial() { return m_pFluid; }
	FESolutesMaterial*	GetSolutesMaterial() { return m_pSolute; }

public: // solutes interface
	int Solutes() override { return m_pSolute->Solutes(); }
	FESolute* GetSolute(int i) override { return m_pSolute->GetSolute(i); }
	vec3d SoluteFlux(FEMaterialPoint& pm, int isol) { return m_pSolute->SoluteFlux(pm, isol); }

private: // material properties
	FEFluid*            m_pFluid;       //!< pointer to fluid material
	FESolutesMaterial*	m_pSolute;      //!< pointer to solutes material

	DECLARE_FECORE_CLASS();
};
