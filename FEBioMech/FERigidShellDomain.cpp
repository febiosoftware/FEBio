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



#include "stdafx.h"
#include "FERigidShellDomain.h"

//-----------------------------------------------------------------------------
FERigidShellDomainOld::FERigidShellDomainOld(FEModel* pfem) : FEElasticShellDomainOld(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything, we need it since 
//       for rigid shell domains we don't want to call the FEElasticShellDomain::Initialize member.
bool FERigidShellDomainOld::Init()
{
	// just call the base class
	return FEShellDomainOld::Init();
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidShellDomainOld::Reset()
{
	// nothing here
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomainOld::StiffnessMatrix(FELinearSystem& LS)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomainOld::InternalForces(FEGlobalVector& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomainOld::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}

//===========================================================================================================

//-----------------------------------------------------------------------------
FERigidShellDomain::FERigidShellDomain(FEModel* pfem) : FEElasticShellDomain(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything, we need it since 
//       for rigid shell domains we don't want to call the FEElasticShellDomain::Initialize member.
bool FERigidShellDomain::Init()
{
	// just call the base class
	return FEShellDomain::Init();
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidShellDomain::Reset()
{
	// nothing here
}

//-----------------------------------------------------------------------------
//! Calculate stiffness contributions for rigid shells.
//! Since rigid elements don't generate stress, we don't need to do
//! anything here.
void FERigidShellDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// Caught you looking!
}


//-----------------------------------------------------------------------------
//! calculate residual forces for rigid shells
//!
void FERigidShellDomain::InternalForces(FEGlobalVector& R)
{
	// Nothing to do.
}

//-----------------------------------------------------------------------------
//! update stresses for rigid shells.
//!
void FERigidShellDomain::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}
