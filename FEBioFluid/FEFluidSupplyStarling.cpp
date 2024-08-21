/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
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



#include "stdafx.h"
#include <sstream>
#include <iostream>
#include <cstdlib>
#include "FEFluidSupplyStarling.h"
#include "FEFluidFSI.h"
//#include "FESolutesMaterialPoint.h"
#include <FECore/FEModel.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFluidSupplyStarling, FEFluidSupply)
	ADD_PARAMETER(m_kp, "kp")->setLongName("filtration coefficient");
	ADD_PARAMETER(m_pv, "pv")->setLongName("external pressure")->setUnits(UNIT_PRESSURE);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEFluidSupplyStarling::FEFluidSupplyStarling(FEModel* pfem) : FEFluidSupply(pfem)
{
	m_kp = 0;
	m_pv = 0;
}

//-----------------------------------------------------------------------------
//! Solvent supply
double FEFluidSupplyStarling::Supply(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& pt = *mp.ExtractData<FEFluidMaterialPoint>();

	// evaluate solvent supply from pressure drop
    double kp = m_kp(mp);
    double pv = m_pv(mp);
	double phiwhat = kp*(pv - pt.m_pf);
/*
	// evaluate solvent supply from concentration drop
	if (mpt) {
		int nsol = mpt->m_nsol;
		for (int isol=0; isol<nsol; ++isol) {
			phiwhat += m_qc[isol]*(m_cv[isol] - mpt->m_c[isol]);
		}
	}
*/
	return phiwhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to strain
mat3d FEFluidSupplyStarling::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return mat3d(mat3dd(0));
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to pressure
double FEFluidSupplyStarling::Tangent_Supply_Dilatation(FEMaterialPoint &mp)
{
    double dpdJ = 0;
    FEFluidFSI* m_pMat = dynamic_cast<FEFluidFSI*>(GetParent());
    if (m_pMat) dpdJ = m_pMat->Fluid()->GetElastic()->Tangent_Strain(mp);
    double kp = m_kp(mp);
	return -kp*dpdJ;
}
/*
//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to concentration
double FEFluidSupplyStarling::Tangent_Supply_Concentration(FEMaterialPoint &mp, const int isol)
{
	FESolutesMaterialPoint& mpt = *mp.ExtractData<FESolutesMaterialPoint>();
	if (isol < mpt.m_nsol) {
		return -m_qc[isol];
	}
	
	return 0;
}
*/
