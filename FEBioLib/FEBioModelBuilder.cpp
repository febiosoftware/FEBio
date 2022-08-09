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
#include "FEBioModelBuilder.h"
#include "FEBioModel.h"
#include <FEBioMech/FEUncoupledMaterial.h>
#include <FEBioMech/FEUDGHexDomain.h>
#include <FEBioMech/FEUT4Domain.h>
#include <FEBioMech/FESSIShellDomain.h>
#include <FEBioMech/RigidBC.h>
#include <FEBioMech/FERigidForce.h>
#include <FEBioMech/FEMechModel.h>

// In FEBio 3, the bulk modulus k must be defined at the top - level.
// However, this could break backward compatibility, so for older file version
// we apply this hack that collects the child moduli and assigns it to the top-level
void FixUncoupledMaterial(FEUncoupledMaterial* mat)
{
	double K = mat->m_K;
	for (int i = 0; i < mat->Properties(); ++i)
	{
		FEUncoupledMaterial* mati = dynamic_cast<FEUncoupledMaterial*>(mat->GetProperty(i));
		if (mati)
		{
			FixUncoupledMaterial(mati);
			K += mati->m_K;
			mati->m_K = 0.0;
		}
	}
	mat->m_K = K;
}

FEBioModelBuilder::FEBioModelBuilder(FEBioModel& fem) : FEModelBuilder(fem)
{

}

void FEBioModelBuilder::AddMaterial(FEMaterial* mat)
{
	FEModel& fem = GetFEModel();
	fem.AddMaterial(mat);

	// For uncoupled materials, we collect the bulk moduli of child materials
	// and assign it to the top-level material (this one)
	FEUncoupledMaterial* pucm = dynamic_cast<FEUncoupledMaterial*>(mat);
	if (pucm) FixUncoupledMaterial(pucm);
}

FEDomain* FEBioModelBuilder::CreateDomain(FE_Element_Spec espec, FEMaterial* mat)
{
	FEModel& fem = GetFEModel();

	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEDomain* pdom = febio.CreateDomain(espec, &fem.GetMesh(), mat);

	// Handle dome special cases
	// TODO: Find a better way of dealing with these special cases
	FEUDGHexDomain* udg = dynamic_cast<FEUDGHexDomain*>(pdom);
	if (udg)
	{
		udg->SetHourGlassParameter(m_udghex_hg);
	}

	FEUT4Domain* ut4 = dynamic_cast<FEUT4Domain*>(pdom);
	if (ut4)
	{
		ut4->SetUT4Parameters(m_ut4_alpha, m_ut4_bdev);
	}

	FESSIShellDomain* ssi = dynamic_cast<FESSIShellDomain*>(pdom);
	if (ssi) {
		ssi->m_bnodalnormals = espec.m_shell_norm_nodal;
	}

	return pdom;
}

//-----------------------------------------------------------------------------
void FEBioModelBuilder::AddRigidComponent(FEStepComponent* pmc)
{
	FEMechModel& fem = static_cast<FEMechModel&>(GetFEModel());

	AddComponent(pmc);

	FERigidFixedBC* prc = dynamic_cast<FERigidFixedBC*>(pmc);
	if (prc) { fem.AddRigidFixedBC(prc); return; }

	FERigidPrescribedBC* prf = dynamic_cast<FERigidPrescribedBC*>(pmc);
	if (prf) { fem.AddRigidPrescribedBC(prf); return; }

	FERigidIC* ric = dynamic_cast<FERigidIC*>(pmc);
	if (ric) { fem.AddRigidInitialCondition(ric); return; }

	FEModelLoad* pml = dynamic_cast<FEModelLoad*>(pmc);
	if (pml) { AddModelLoad(pml); return; }

	assert(false);
}
