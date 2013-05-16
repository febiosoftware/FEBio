#include "stdafx.h"
#include "FEBioContactSection.h"
#include <FEBioLib/FESlidingInterface.h>
#include <FEBioLib/FESlidingInterface2.h>
#include <FEBioLib/FESlidingInterface3.h>
#include <FEBioLib/FETiedInterface.h>
#include <FEBioLib/FETiedBiphasicInterface.h>
#include <FEBioLib/FEFacet2FacetSliding.h>
#include <FEBioLib/FESlidingInterfaceBW.h>
#include <FEBioLib/FEPeriodicBoundary.h>
#include <FEBioLib/FESurfaceConstraint.h>
#include <FEBioLib/FERigidWallInterface.h>
#include <FEBioLib/FERigidJoint.h>
#include <FEBioLib/FEAugLagLinearConstraint.h>

//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.x
	int nversion = m_pim->Version();
	if (nversion < 0x0200) throw XMLReader::InvalidTag(tag);

	// loop over tags
	++tag;
	do
	{
		if (tag == "contact")
		{
			// get the contact type
			const char* sztype = tag.AttributeValue("type");
			if      (strcmp(sztype, "sliding_with_gaps"     ) == 0) ParseSlidingInterface     (tag);
			else if (strcmp(sztype, "facet-to-facet sliding") == 0) ParseFacetSlidingInterface(tag);
			else if (strcmp(sztype, "sliding2"              ) == 0) ParseSlidingInterface2    (tag);
			else if (strcmp(sztype, "sliding3"              ) == 0) ParseSlidingInterface3    (tag);
			else if (strcmp(sztype, "tied"                  ) == 0) ParseTiedInterface        (tag);
			else if (strcmp(sztype, "periodic boundary"     ) == 0) ParsePeriodicBoundary     (tag);
			else if (strcmp(sztype, "surface constraint"    ) == 0) ParseSurfaceConstraint    (tag);
			else if (strcmp(sztype, "rigid_wall"            ) == 0) ParseRigidWall            (tag);
			else if (strcmp(sztype, "rigid"                 ) == 0) ParseRigidInterface       (tag);
			else if (strcmp(sztype, "rigid joint"           ) == 0) ParseRigidJoint           (tag);
			else if (strcmp(sztype, "linear constraint"     ) == 0) ParseLinearConstraint     (tag);
			else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- S L I D I N G   W I T H   G A P S ---
void FEBioContactSection::ParseSlidingInterface(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FESlidingInterface* ps = new FESlidingInterface(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		// read parameters
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FESlidingSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, true);
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- F A C E T   T O   F A C E T   S L I D I N G ---
void FEBioContactSection::ParseFacetSlidingInterface(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FEFacet2FacetSliding* ps = new FEFacet2FacetSliding(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		// read parameters
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FEFacetSlidingSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, false);
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- S L I D I N G   I N T E R F A C E   2 ---
void FEBioContactSection::ParseSlidingInterface2(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FESlidingInterface2* ps = new FESlidingInterface2(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		// read parameters
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FESlidingSurface2& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, false);
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- S L I D I N G   I N T E R F A C E   3 ---
void FEBioContactSection::ParseSlidingInterface3(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FESlidingInterface3* ps = new FESlidingInterface3(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();
	
	++tag;
	do
	{
		// read parameters
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;
				
				FESlidingSurface3& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);
				
				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}
				
				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, false);
			}
			else throw XMLReader::InvalidTag(tag);
		}
		
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- T I E D   C O N T A C T  ---
void FEBioContactSection::ParseTiedInterface(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FETiedInterface* ps = new FETiedInterface(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FETiedContactSurface& s = (ntype == 1? ps->ms : ps->ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, true);
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- P E R I O D I C   B O U N D A R Y  ---
void FEBioContactSection::ParsePeriodicBoundary(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FEPeriodicBoundary* ps = new FEPeriodicBoundary(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FEPeriodicSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, true);
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- S U R F A C E   C O N S T R A I N T ---
void FEBioContactSection::ParseSurfaceConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FESurfaceConstraint* ps = new FESurfaceConstraint(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FESurfaceConstraintSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);
				m.AddSurface(&s);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, true);
			}
			else throw XMLReader::InvalidTag(tag);
		}

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- R I G I D   W A L L   I N T E R F A C E ---
void FEBioContactSection::ParseRigidWall(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FERigidWallInterface* ps = new FERigidWallInterface(&fem);
	fem.AddContactInterface(ps);

	FEParameterList& pl = ps->GetParameterList();

	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			if (tag == "plane")
			{
				ps->SetMasterSurface(new FEPlane(&fem));
				FEPlane& pl = dynamic_cast<FEPlane&>(*ps->m_mp);
				const char* sz = tag.AttributeValue("lc", true);
				if (sz)	pl.m_nplc = atoi(sz) - 1;

				double* a = pl.GetEquation();
				tag.value(a, 4);
			}
			else if (tag == "sphere")
			{
				ps->SetMasterSurface(new FERigidSphere(&fem));
				FERigidSphere& s = dynamic_cast<FERigidSphere&>(*ps->m_mp);
				++tag;
				do
				{
					if      (tag == "center") tag.value(s.m_rc);
					else if (tag == "radius") tag.value(s.m_R);
					else if (tag == "xtrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[0] = atoi(szlc) - 1;
					}
					else if (tag == "ytrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[1] = atoi(szlc) - 1;
					}
					else if (tag == "ztrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[2] = atoi(szlc) - 1;
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}
			else if (tag == "surface")
			{
				FERigidWallSurface& s = ps->m_ss;

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt, true);
			}
			else throw XMLReader::InvalidTag(tag);
		}
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
// --- R I G I D   B O D Y   I N T E R F A C E ---
void FEBioContactSection::ParseRigidInterface(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// count how many rigid nodes there are
	int nrn= 0;
	XMLTag t(tag); ++t;
	while (!t.isend()) { nrn++; ++t; }

	++tag;
	int id, rb;
	for (int i=0; i<nrn; ++i)
	{
		id = atoi(tag.AttributeValue("id"))-1;
		rb = atoi(tag.AttributeValue("rb"))-1;

		FERigidNode* prn = new FERigidNode;

		prn->nid = id;
		prn->rid = rb;
		fem.AddRigidNode(prn);

		if (m_pim->m_nsteps > 0)
		{
			GetStep()->AddBoundaryCondition(prn);
			prn->Deactivate();
		}

		++tag;
	}
}

//-----------------------------------------------------------------------------
// --- R I G I D   J O I N T   I N T E R F A C E ---
void FEBioContactSection::ParseRigidJoint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	FERigidJoint* prj = new FERigidJoint(&fem);
	FEParameterList& pl = prj->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
	prj->m_nRBa--;
	prj->m_nRBb--;
	fem.AddNonlinearConstraint(prj);
}

//-----------------------------------------------------------------------------
// --- L I N E A R   C O N S T R A I N T ---
void FEBioContactSection::ParseLinearConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	// create a new linear constraint manager
	FELinearConstraintSet* pLCS = new FELinearConstraintSet(&fem);
	fem.AddNonlinearConstraint(pLCS);

	// read the linear constraints
	++tag;
	do
	{
		if (tag == "linear_constraint")
		{
			FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;

			FEAugLagLinearConstraint::DOF dof;
			++tag;
			do
			{
				if (tag == "node")
				{
					tag.value(dof.val);
					int node;
					tag.AttributeValue("id", node);
					dof.node = node - 1;

					const char* szbc = tag.AttributeValue("bc");
					if      (strcmp(szbc, "x") == 0) dof.bc = 0;
					else if (strcmp(szbc, "y") == 0) dof.bc = 1;
					else if (strcmp(szbc, "z") == 0) dof.bc = 2;
					else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

					pLC->m_dof.push_back(dof);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			// add the linear constraint to the system
			pLCS->add(pLC);
		}
		else if (tag == "tol"    ) tag.value(pLCS->m_tol);
		else if (tag == "penalty") tag.value(pLCS->m_eps);
		else if (tag == "maxaug") tag.value(pLCS->m_naugmax);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioContactSection::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt, bool bnodal)
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
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
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
				int ne[4];
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
