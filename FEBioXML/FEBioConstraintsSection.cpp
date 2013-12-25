#include "stdafx.h"
#include "FEBioConstraintsSection.h"
#include "FEBioMech/FERigidMaterial.h"
#include "FEBioMech/FEPointConstraint.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"

//=============================================================================
//
//                  C O N S T R A I N T S   S E C T I O N
//
//=============================================================================

void FEBioConstraintsSection::Parse(XMLTag &tag)
{
	// make sure there is something to read
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "rigid_body") ParseRigidConstraint(tag);
		else if (tag == "constraint")
		{
			const char* sztype = tag.AttributeValue("type");
			FENLConstraint* plc = fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, sztype, m_pim->GetFEModel());
			if (plc == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			FEParameterList& pl = plc->GetParameterList();

			++tag;
			do
			{
				if (m_pim->ReadParameter(tag, pl) == false)
				{
					if (tag == "surface")
					{
						const char* sztype = tag.AttributeValue("type");
						FESurface* ps = plc->GetSurface(sztype);
						if (ps == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
						ParseSurfaceSection(tag, *ps, 0, true);
					}
					else throw XMLReader::InvalidTag(tag);
				}
				++tag;
			}
			while (!tag.isend());

			FEModel& fem = *GetFEModel();
			fem.AddNonlinearConstraint(plc);

			// add this boundary condition to the current step
			if (m_pim->m_nsteps > 0)
			{
				GetStep()->AddConstraint(plc);
				plc->Deactivate();
			}
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
	FEAnalysis* pStep = (m_pim->m_nsteps > 0 ? GetStep() : 0);

	const char* szm = tag.AttributeValue("mat");
	assert(szm);

	// get the material ID
	int nmat = atoi(szm);
	if ((nmat <= 0) || (nmat > fem.Materials())) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

	// make sure this is a valid rigid material
	FEMaterial* pm = fem.GetMaterial(nmat-1);
	if (pm->IsRigid() == false) throw XMLReader::InvalidAttributeValue(tag, "mat", szm);

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

				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement(&fem);
				pDC->id = nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);

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

				FERigidBodyForce* pFC = new FERigidBodyForce(&fem);
				pFC->id = nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else if (strcmp(szt, "fixed") == 0)
			{
				FERigidBodyFixedBC* pBC = new FERigidBodyFixedBC(&fem);
				pBC->id = nmat;
				pBC->bc = bc;
				fem.m_RBC.push_back(pBC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RBC.size()-1;
					FERigidBodyFixedBC* pBC = fem.m_RBC[n];
					pStep->AddBoundaryCondition(pBC);
					pBC->Deactivate();
				}
			}
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

				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement(&fem);
				pDC->id = nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);

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

				FERigidBodyForce* pFC = new FERigidBodyForce(&fem);
				pFC->id = nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else if (strcmp(szt, "fixed") == 0)
			{
				FERigidBodyFixedBC* pBC = new FERigidBodyFixedBC(&fem);
				pBC->id = nmat;
				pBC->bc = bc;
				fem.m_RBC.push_back(pBC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RBC.size()-1;
					FERigidBodyFixedBC* pBC = fem.m_RBC[n];
					pStep->AddBoundaryCondition(pBC);
					pBC->Deactivate();
				}
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioConstraintsSection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	// count nr of faces
	int faces = 0, N, nf[8];
	XMLTag t(tag); ++t;
	while (!t.isend()) { faces++; ++t; }

	// allocate storage for faces
	s.create(faces);

	// read faces
	++tag;
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		// set the element type/integration rule
		if (bnodal)
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4NI);
			else if (tag == "tri3" ) el.SetType(FE_TRI3NI );
			else if (tag == "tri6" ) el.SetType(FE_TRI6NI);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else throw XMLReader::InvalidTag(tag);
		}

		N = el.Nodes();

		if (nfmt == 0)
		{
			tag.value(nf, N);
			for (int j=0; j<N; ++j) 
			{
				int nid = nf[j]-1;
				if ((nid<0)||(nid>= NN)) throw XMLReader::InvalidValue(tag);
				el.m_node[j] = nid;
			}
		}
		else if (nfmt == 1)
		{
			tag.value(nf, 2);
			FEElement* pe = m.FindElementFromID(nf[0]);
			if (pe)
			{
				int ne[8];
				int nn = m.GetFace(*pe, nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
				el.m_nelem = nf[0];
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}
