#include "stdafx.h"
#include "FEBioContactSection.h"
#include "FEBioMech/FESlidingInterface.h"
#include "FEBioMix/FESlidingInterface2.h"
#include "FEBioMix/FESlidingInterface3.h"
#include "FEBioMix/FETiedBiphasicInterface.h"
#include "FEBioMech/FETiedInterface.h"
#include "FEBioMech/FEFacet2FacetTied.h"
#include "FEBioMech/FEFacet2FacetSliding.h"
#include "FEBioMech/FESlidingInterfaceBW.h"
#include "FEBioMech/FEPeriodicBoundary.h"
#include "FEBioMech/FESurfaceConstraint.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FERigidJoint.h"
#include "FEBioMech/FEAugLagLinearConstraint.h"

//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.x
	int nversion = m_pim->Version();
	if (nversion < 0x0200) throw XMLReader::InvalidTag(tag);

	FEModel& fem = *GetFEModel();

	// loop over tags
	++tag;
	do
	{
		if (tag == "contact")
		{
			// get the contact type
			const char* sztype = tag.AttributeValue("type");

			// Not all contact interfaces can be automated, so we first handle these special cases
			if      (strcmp(sztype, "rigid_wall"            ) == 0) ParseRigidWall            (tag);
			else if (strcmp(sztype, "rigid"                 ) == 0) ParseRigidInterface       (tag);
			else if (strcmp(sztype, "rigid joint"           ) == 0) ParseRigidJoint           (tag);
			else if (strcmp(sztype, "linear constraint"     ) == 0) ParseLinearConstraint     (tag);
			else 
			{
				// If we get here, we try to create a contact interface
				// using the FEBio kernel. 
				FEBioKernel& febio = FEBioKernel::GetInstance();
				FEContactInterface* pci = febio.Create<FEContactInterface>(sztype, GetFEModel());
				if (pci)
				{
					fem.AddSurfacePairInteraction(pci);
					ParseContactInterface(tag, pci);
				}
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioContactSection::ParseContactInterface(XMLTag& tag, FEContactInterface* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();

	// read the parameters
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

				FESurface& s = *(ntype == 1? pci->GetMasterSurface() : pci->GetSlaveSurface());
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
	fem.AddSurfacePairInteraction(ps);

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

	int N, nf[8];

	// count nr of faces
	int faces = tag.children();

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
