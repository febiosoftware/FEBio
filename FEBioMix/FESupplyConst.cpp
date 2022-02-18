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
#include "FESupplyConst.h"

// define the material parameters
BEGIN_FECORE_CLASS(FESupplyConst, FESoluteSupply)
	ADD_PARAMETER(m_supp, "supp");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FESupplyConst::FESupplyConst(FEModel* pfem) : FESoluteSupply(pfem)
{
	m_supp = 0;
}

//-----------------------------------------------------------------------------
//! Solute supply
double FESupplyConst::Supply(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to strain
double FESupplyConst::Tangent_Supply_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of solute supply with respect to concentration
double FESupplyConst::Tangent_Supply_Concentration(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand complex supply
double FESupplyConst::ReceptorLigandSupply(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Solute supply at steady-state
double FESupplyConst::SupplySS(FEMaterialPoint& mp)
{
	// --- constant solubility ---
	
	return m_supp;
}

//-----------------------------------------------------------------------------
//! Receptor-ligand concentration at steady-state
double FESupplyConst::ReceptorLigandConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Referential solid supply (moles of solid/referential volume/time)
double FESupplyConst::SolidSupply(FEMaterialPoint& mp)
{
	return ReceptorLigandSupply(mp);
}

//-----------------------------------------------------------------------------
//! Referential solid concentration (moles of solid/referential volume)
//! at steady-state
double FESupplyConst::SolidConcentrationSS(FEMaterialPoint& mp)
{
	return 0;
}


