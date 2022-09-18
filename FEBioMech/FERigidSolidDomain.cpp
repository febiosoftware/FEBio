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
#include "FERigidSolidDomain.h"
#include "FERigidMaterial.h"
#include <FECore/FEMesh.h>

//-----------------------------------------------------------------------------
FERigidSolidDomain::FERigidSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) {}

//-----------------------------------------------------------------------------
// NOTE: Although this function doesn't do anything we need it because we don't
//       want to call the FEElasticSolidDomain::Initialize function.
bool FERigidSolidDomain::Init()
{
	if (FESolidDomain::Init() == false) return false;

	// set the rigid nodes ID
	// Altough this is already done in FERigidSystem::CreateObjects()
	// we do it here again since the mesh adaptors will call this function
	// and they may have added some nodes
	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(m_pMat); assert(pm);
	if (pm == nullptr) return false;
	int rbid = pm->GetRigidBodyID();
	for (int i = 0; i < Nodes(); ++i)
	{
		FENode& node = Node(i);
		node.m_rid = rbid;
	}
	return true;
}

//-----------------------------------------------------------------------------
// We need to override it since the base class version will not work for rigid domains.
void FERigidSolidDomain::Reset()
{
	// nothing to reset here
}

//-----------------------------------------------------------------------------
//! Calculates the stiffness matrix for 3D rigid elements.
//! Rigid elements don't generate stress, so there is nothing to do here
void FERigidSolidDomain::StiffnessMatrix(FELinearSystem& LS)
{
	// Caught you looking!
}

//-----------------------------------------------------------------------------
// Rigid bodies do not generate stress so there is nothing to do here
void FERigidSolidDomain::InternalForces(FEGlobalVector& R)
{
	// what you looking at ?!
}

//-----------------------------------------------------------------------------
//! We don't need to update the stresses for rigid elements
//!
void FERigidSolidDomain::Update(const FETimeInfo& tp)
{
	// Nothing to see here. Please move on.
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::MassMatrix(FELinearSystem& LS, double scale)
{
	// Only crickets here ... 
}

//-----------------------------------------------------------------------------
void FERigidSolidDomain::InertialForces(FEGlobalVector& R, std::vector<double>& F)
{
	// chirp, chirp ...
}
