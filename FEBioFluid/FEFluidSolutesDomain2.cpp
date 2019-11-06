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
#include "FEFluidSolutesDomain2.h"
#include "FECore/log.h"
#include "FECore/DOFS.h"
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/sys.h>
#include "FEBioFluidSolutes.h"
#include "FEFluidSolutesMaterial2.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
//! constructor
//! Some derived classes will pass 0 to the pmat, since the pmat variable will be
//! to initialize another material. These derived classes will set the m_pMat variable as well.
FEFluidSolutesDomain2::FEFluidSolutesDomain2(FEModel* pfem) : FESolidDomain(pfem), FEFluidDomain3D(pfem), FESolutesDomain(pfem)
{
	m_pMat = 0;
	m_activeDomain = 0;
}


//-----------------------------------------------------------------------------
FEFluidSolutesDomain2::~FEFluidSolutesDomain2()
{

}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::SetActiveDomain(int n)
{
	m_activeDomain = n;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::Serialize(DumpStream& ar)
{
	ar & m_activeDomain;
	if (m_activeDomain == FLUID_DOMAIN)	FEFluidDomain3D::Serialize(ar);
	else FESolutesDomain::Serialize(ar);
}

//-----------------------------------------------------------------------------
const FEDofList& FEFluidSolutesDomain2::GetDOFList() const
{
	if (m_activeDomain == FLUID_DOMAIN) return FEFluidDomain3D::GetDOFList();
	else return FESolutesDomain::GetDOFList();
}

//-----------------------------------------------------------------------------
//! Assign material
void FEFluidSolutesDomain2::SetMaterial(FEMaterial* pmat)
{
	m_pMat = dynamic_cast<FEFluidSolutesMaterial2*>(pmat);
	assert(m_pMat);

	// set the fluid material first
	FEFluidDomain3D::SetMaterial(m_pMat->GetFluidMaterial());
	FESolutesDomain::SetMaterial(m_pMat->GetSolutesMaterial());

	// Since the base class is shared, we need to call the base class last
	// to make sure that the base class knows this material is the correct one.
	FEDomain::SetMaterial(pmat);
}

//-----------------------------------------------------------------------------
bool FEFluidSolutesDomain2::Init()
{
	if (FEFluidDomain3D::Init() == false) return false;
	if (FESolutesDomain::Init() == false) return false;
	return true;
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::Reset()
{
	// reset base class
	FEFluidDomain3D::Reset();
	FESolutesDomain::Reset();
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::Activate()
{
	SetActiveDomain(FLUID_DOMAIN);
	FEFluidDomain3D::Activate();
	FESolutesDomain::Activate();
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::InitMaterialPoints()
{
	FEFluidDomain3D::InitMaterialPoints();
	FESolutesDomain::InitMaterialPoints();
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::PreSolveUpdate(const FETimeInfo& tp) 
{ 
	if (m_activeDomain == 0) FEFluidDomain3D::PreSolveUpdate(tp);
	else FESolutesDomain::PreSolveUpdate(tp);
}

//-----------------------------------------------------------------------------
void FEFluidSolutesDomain2::Update(const FETimeInfo& tp)
{ 
	if (m_activeDomain == 0) FEFluidDomain3D::Update(tp);
	else FESolutesDomain::Update(tp);
}
