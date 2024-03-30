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
#include "FERigidRotationVector.h"

BEGIN_FECORE_CLASS(FERigidRotationVector, FERigidBC)
	ADD_PARAMETER(m_rigidMat, "rb")->setEnums("$(rigid_materials)")->setLongName("Rigid material");
	ADD_PARAMETER(m_vx, "vx")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_RADIAN);
	ADD_PARAMETER(m_vy, "vy")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_RADIAN);
	ADD_PARAMETER(m_vz, "vz")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE)->setUnits(UNIT_RADIAN);
END_FECORE_CLASS();

FERigidRotationVector::FERigidRotationVector(FEModel* fem) : FERigidBC(fem)
{
	m_vx = m_vy = m_vz = 0.0;

	m_rc[0] = new FERigidPrescribedBC(fem); m_rc[0]->SetBC(3);
	m_rc[1] = new FERigidPrescribedBC(fem); m_rc[1]->SetBC(4);
	m_rc[2] = new FERigidPrescribedBC(fem); m_rc[2]->SetBC(5);
}

FERigidRotationVector::~FERigidRotationVector()
{
	delete m_rc[0];
	delete m_rc[1];
	delete m_rc[2];
}

bool FERigidRotationVector::Init()
{
	for (int i = 0; i < 3; ++i)
	{
		m_rc[i]->SetRigidMaterial(m_rigidMat);
		if (m_rc[i]->Init() == false) return false;
	}

	return true;
}

void FERigidRotationVector::Activate()
{
	FERigidBC::Activate();
	m_rc[0]->Activate();
	m_rc[1]->Activate();
	m_rc[2]->Activate();

	m_rc[0]->SetValue(m_vx);
	m_rc[1]->SetValue(m_vy);
	m_rc[2]->SetValue(m_vz);
}

void FERigidRotationVector::Deactivate()
{
	FERigidBC::Deactivate();
	m_rc[0]->Deactivate();
	m_rc[1]->Deactivate();
	m_rc[2]->Deactivate();
}

void FERigidRotationVector::InitTimeStep()
{
	m_rc[0]->SetValue(m_vx);
	m_rc[1]->SetValue(m_vy);
	m_rc[2]->SetValue(m_vz);
}

void FERigidRotationVector::Serialize(DumpStream& ar)
{
	FERigidBC::Serialize(ar);

	if (ar.IsShallow() == false)
	{
		for (int j = 0; j < 3; ++j)
		{
			// We need to serialize a reference to avoid that
			// restart will try to allocate the m_rc variables
			FERigidPrescribedBC& rc = *m_rc[j];
			ar & rc;
		}
	}
}
