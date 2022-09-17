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
#include "FERigidMaterial.h"
#include "FERigidBody.h"
#include "FEMechModel.h"
#include <FECore/log.h>

// define the material parameters
BEGIN_FECORE_CLASS(FERigidMaterial, FESolidMaterial)
	ADD_PARAMETER(m_E      , FE_RANGE_GREATER_OR_EQUAL(0.0), "E"         )->setLongName("Young's modulus");
	ADD_PARAMETER(m_v      , FE_RANGE_RIGHT_OPEN(-1.0, 0.5), "v")->setLongName("Poisson's ratio");
	ADD_PARAMETER(m_com    , "override_com")->SetFlags(FE_PARAM_WATCH);
	ADD_PARAMETER(m_rc    , "center_of_mass")->SetWatchVariable(&m_com);
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// constructor
FERigidMaterial::FERigidMaterial(FEModel* pfem) : FESolidMaterial(pfem)
{
	m_E = 1;
	m_v = 0;
	m_rc = vec3d(0, 0, 0);

	m_com = false;	// calculate COM automatically
	m_pmid = -1;
	m_nRB = -1;

	m_binit = false;
}

//-----------------------------------------------------------------------------
// Initialize rigid material data
bool FERigidMaterial::Init()
{
	if (FESolidMaterial::Init() == false) return false;

	if (m_binit == false)
	{
		// get this rigid body's ID
		FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
		FERigidBody& rb = *fem.GetRigidBody(GetRigidBodyID());

		// only set the rigid body com if this is the main rigid body material
		if (rb.GetMaterialID() == GetID()-1)
		{
			// first, calculate the mass
			rb.UpdateMass();

			// next, calculate the center of mass, or just set it
			if (m_com == false)
			{
				rb.UpdateCOM();
			}
			else
			{
				rb.SetCOM(m_rc); 
			}

			// finally, determine moi
			rb.UpdateMOI();
		}

		if (m_pmid  > -1)
		{
			FERigidMaterial* ppm = dynamic_cast<FERigidMaterial*>(GetFEModel()->GetMaterial(m_pmid-1));
			if (ppm == 0) {
				feLogError("parent of rigid material %s is not a rigid material\n", GetName().c_str());
				return false;
			}

			FERigidBody& prb = *fem.GetRigidBody(ppm->GetRigidBodyID());
			rb.m_prb = &prb;

			// mark all degrees of freedom as prescribed
			rb.m_BC[0] = DOF_PRESCRIBED;
			rb.m_BC[1] = DOF_PRESCRIBED;
			rb.m_BC[2] = DOF_PRESCRIBED;
			rb.m_BC[3] = DOF_PRESCRIBED;
			rb.m_BC[4] = DOF_PRESCRIBED;
			rb.m_BC[5] = DOF_PRESCRIBED;
		}

		m_binit = true;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Serialize data to or from the dump file
void FERigidMaterial::Serialize(DumpStream &ar)
{
	FESolidMaterial::Serialize(ar);
	ar & m_com;
	ar & m_nRB;
}
