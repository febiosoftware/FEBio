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
#include "stdafx.h"
#include "FEFluidSolutesMaterial2.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FEFluidSolutesMaterial2, FEMaterial)
	ADD_PROPERTY(m_pFluid, "fluid");
	ADD_PROPERTY(m_pSolute, "solutes");
END_FECORE_CLASS();

FEFluidSolutesMaterial2::FEFluidSolutesMaterial2(FEModel* fem) : FEMaterial(fem)
{
	m_pFluid = nullptr;
	m_pSolute = nullptr;
}

// returns a pointer to a new material point object
FEMaterialPoint* FEFluidSolutesMaterial2::CreateMaterialPointData()
{
	FEMaterialPoint* fp = m_pFluid->CreateMaterialPointData();
	FEMaterialPoint* sp = m_pSolute->CreateMaterialPointData();
	fp->SetNext(sp);
	return fp;
}
