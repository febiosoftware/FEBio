#include "stdafx.h"
#include "FEBioConstraintsSection.h"
#include <FECore/FERigid.h>
#include <FEBioLib/FEPointConstraint.h>

//=============================================================================
//
//                  C O N S T R A I N T S   S E C T I O N
//
//=============================================================================

void FEBioConstraintsSection::Parse(XMLTag &tag)
{
	// This section is only allowed in the new format
	if (m_pim->Version() < 0x0101) throw XMLReader::InvalidTag(tag);

	// make sure there is something to read
	if (tag.isleaf()) return;

	// get the FEBio kernel
	FEBioKernel& febio = FEBioKernel::GetInstance();

	++tag;
	do
	{
		if (tag == "rigid_body") ParseRigidConstraint(tag);
		else if (tag == "point") ParsePointConstraint(tag);
		else if (tag == "constraint")
		{
			const char* sztype = tag.AttributeValue("type");
			FENLConstraint* plc = febio.Create<FENLConstraint>(sztype, m_pim->GetFEModel());
			if (plc == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			FEParameterList& pl = plc->GetParameterList();

			++tag;
			do
			{
				if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			FEModel& fem = *GetFEModel();
			fem.AddNonlinearConstraint(plc);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioConstraintsSection::ParseRigidConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEAnalysisStep* pStep = (m_pim->m_nsteps > 0 ? GetStep() : 0);

	const char* szm = tag.AttributeValue("mat");
	assert(szm);

	int nmat = atoi(szm);
	if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(nmat-1));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	++tag;
	do
	{
		if (strncmp(tag.Name(), "trans_", 6) == 0)
		{
			const char* szt = tag.AttributeValue("type");

			int bc = -1;
			if      (tag.Name()[6] == 'x') bc = 0;
			else if (tag.Name()[6] == 'y') bc = 1;
			else if (tag.Name()[6] == 'z') bc = 2;
			assert(bc >= 0);
			
			if (strcmp(szt, "prescribed") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc) - 1;

				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement;
				pDC->id = nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);
				pm->m_bc[bc] = lc+1;

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RDC.size()-1;
					FERigidBodyDisplacement* pDC = fem.m_RDC[n];
					pStep->AddBoundaryCondition(pDC);
					pDC->Deactivate();
				}
			}
			else if (strcmp(szt, "force") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc) - 1;

				FERigidBodyForce* pFC = new FERigidBodyForce;
				pFC->id = nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);
				pm->m_bc[bc] = 0;

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else if (strcmp(szt, "fixed") == 0) pm->m_bc[bc] = -1;
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
		else if (strncmp(tag.Name(), "rot_", 4) == 0)
		{
			const char* szt = tag.AttributeValue("type");

			int bc = -1;
			if      (tag.Name()[4] == 'x') bc = 3;
			else if (tag.Name()[4] == 'y') bc = 4;
			else if (tag.Name()[4] == 'z') bc = 5;
			assert(bc >= 0);

			if (strcmp(szt, "prescribed") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc) - 1;

				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement;
				pDC->id = nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);
				pm->m_bc[bc] = lc+1;

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RDC.size()-1;
					FERigidBodyDisplacement* pDC = fem.m_RDC[n];
					pStep->AddBoundaryCondition(pDC);
					pDC->Deactivate();
				}
			}
			else if (strcmp(szt, "force") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc) - 1;

				FERigidBodyForce* pFC = new FERigidBodyForce;
				pFC->id = nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);
				pm->m_bc[bc] = 0;

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else if (strcmp(szt, "fixed") == 0) pm->m_bc[bc] = -1;
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioConstraintsSection::ParsePointConstraint(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	int node = -1;
	double	eps;

	++tag;
	do
	{
		if (tag == "node") 
		{
			tag.value(node);
			if (node <= 0) throw XMLReader::InvalidValue(tag);
		}
		else if (tag == "penalty") tag.value(eps);
		++tag;
	}
	while (!tag.isend());
	if (node == -1) throw XMLReader::Error();

	FEPointConstraint* pc = new FEPointConstraint(&fem);
	pc->m_eps = eps;
	pc->m_node = node-1;
	fem.AddNonlinearConstraint(pc);
}
