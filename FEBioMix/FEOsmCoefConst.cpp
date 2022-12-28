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
#include "FEOsmCoefConst.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_FECORE_CLASS(FEOsmCoefConst, FEOsmoticCoefficient)
	ADD_PARAMETER(m_osmcoef, FE_RANGE_GREATER_OR_EQUAL(0.0), "osmcoef")->setLongName("osmotic coefficient");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! Constructor. 
FEOsmCoefConst::FEOsmCoefConst(FEModel* pfem) : FEOsmoticCoefficient(pfem)
{
	m_osmcoef = 1;
}

//-----------------------------------------------------------------------------
//! Osmotic coefficient
double FEOsmCoefConst::OsmoticCoefficient(FEMaterialPoint& mp)
{
	// --- constant osmotic coefficient ---
	
	return m_osmcoef;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to strain
double FEOsmCoefConst::Tangent_OsmoticCoefficient_Strain(FEMaterialPoint &mp)
{
	return 0;
}

//-----------------------------------------------------------------------------
//! Tangent of osmotic coefficient with respect to concentration
double FEOsmCoefConst::Tangent_OsmoticCoefficient_Concentration(FEMaterialPoint &mp, const int isol)
{
	return 0;
}

