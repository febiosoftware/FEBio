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
#include "FESolventSupplyStarling.h"
#include "FESolutesMaterialPoint.h"
#include "FECore/FEModel.h"

// define the material parameters
BEGIN_FECORE_CLASS(FESolventSupplyStarling, FESolventSupply)
	ADD_PARAMETER(m_kp, "kp")->setLongName("filtration coefficient");
	ADD_PARAMETER(m_pv, "pv")->setLongName("external pressure");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESolventSupplyStarling::FESolventSupplyStarling(FEModel* pfem) : FESolventSupply(pfem)
{
	m_kp = 0;
	m_pv = 0;

    // get number of DOFS
	DOFS& fedofs = pfem->GetDOFS();
    int MAX_CDOFS = fedofs.GetVariableSize("concentration");
    
    if (MAX_CDOFS > 0) {
        m_qc.assign(MAX_CDOFS,0);
        m_cv.assign(MAX_CDOFS,0);
    }
}

//-----------------------------------------------------------------------------
//! Solvent supply
double FESolventSupplyStarling::Supply(FEMaterialPoint& mp)
{
	FEBiphasicMaterialPoint& ppt = *mp.ExtractData<FEBiphasicMaterialPoint>();
	FESolutesMaterialPoint* mpt = mp.ExtractData<FESolutesMaterialPoint>();

	// evaluate solvent supply from pressure drop
	double phiwhat = m_kp*(m_pv - ppt.m_p);
	
	// evaluate solvent supply from concentration drop
	if (mpt) {
		int nsol = mpt->m_nsol;
		for (int isol=0; isol<nsol; ++isol) {
			phiwhat += m_qc[isol]*(m_cv[isol] - mpt->m_c[isol]);
		}
	}
	
	return phiwhat;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to strain
mat3ds FESolventSupplyStarling::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	mat3dd Phie(Supply(mp));
	
	return Phie;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to pressure
double FESolventSupplyStarling::Tangent_Supply_Pressure(FEMaterialPoint &mp)
{
	return -m_kp;
}

//-----------------------------------------------------------------------------
//! Tangent of solvent supply with respect to concentration
double FESolventSupplyStarling::Tangent_Supply_Concentration(FEMaterialPoint &mp, const int isol)
{
	FESolutesMaterialPoint& mpt = *mp.ExtractData<FESolutesMaterialPoint>();
	if (isol < mpt.m_nsol) {
		return -m_qc[isol];
	}
	
	return 0;
}

