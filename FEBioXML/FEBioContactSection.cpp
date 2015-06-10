#include "stdafx.h"
#include "FEBioContactSection.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEAugLagLinearConstraint.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection::Parse(XMLTag& tag)
{
	// make sure that the version is 2.x
	int nversion = m_pim->Version();
	if (nversion < 0x0200) throw XMLReader::InvalidTag(tag);

	// make sure there are children
	if (tag.isleaf()) return;

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
			if      (strcmp(sztype, "rigid_wall"             ) == 0) ParseRigidWall            (tag);
			else if (strcmp(sztype, "rigid"                  ) == 0) ParseRigidInterface       (tag);
			else if (strcmp(sztype, "linear constraint"      ) == 0) ParseLinearConstraint     (tag);
			else 
			{
				// If we get here, we try to create a contact interface
				// using the FEBio kernel. 
				FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fecore_new<FESurfacePairInteraction>(FESURFACEPAIRINTERACTION_ID, sztype, &fem));
				if (pci)
				{
					fem.AddSurfacePairInteraction(pci);
					ParseContactInterface(tag, pci);
					// add this boundary condition to the current step
					if (m_pim->m_nsteps > 0)
					{
						GetStep()->AddSurfacePairInteraction(pci);
						pci->Deactivate();
					}
				}
				else
				{
					FEModel& fem = *GetFEModel();

					// Some constraints were initially defined in the Contact section, although
					// now it is preferred that they are defined in the Constraints section. For backward
					// compatibility we still allow constraints to be defined in this section. 
					FENLConstraint* pc = fecore_new<FENLConstraint>(FENLCONSTRAINT_ID, sztype, &fem);
					if (pc)
					{
						FEParameterList& pl = pc->GetParameterList();
						++tag;
						do
						{
							if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
							++tag;
						}
						while (!tag.isend());
						fem.AddNonlinearConstraint(pc);
						if (m_pim->m_nsteps > 0)
						{
							GetStep()->AddConstraint(pc);
							pc->Deactivate();
						}
					}
					else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
				}
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioContactSection::ParseContactInterface(XMLTag& tag, FESurfacePairInteraction* pci)
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
				int ntype = 0;
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

				// see if the set attribute is defined
				const char* szset = tag.AttributeValue("set", true);
				if (szset)
				{
					// make sure this tag does not have any children
					if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);

					// see if we can find the facet set
					FEFacetSet* ps = 0;
					for (int i=0; i<m.FacetSets(); ++i)
					{
						FEFacetSet& fi = m.FacetSet(i);
						if (strcmp(fi.GetName(), szset) == 0)
						{
							ps = &fi;
							break;
						}
					}

					// create a surface from the facet set
					if (ps)
					{
						if (BuildSurface(s, *ps, pci->UseNodalIntegration()) == false) throw XMLReader::InvalidTag(tag);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
				}
				else 
				{
					// read the surface section
					if (ParseSurfaceSection(tag, s, nfmt, pci->UseNodalIntegration()) == false) throw XMLReader::InvalidTag(tag);
				}
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

	FERigidWallInterface* ps = new FERigidWallInterface(&fem);
	fem.AddSurfacePairInteraction(ps);

	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, ps) == false)
		{
			if (tag == "surface")
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

	int NMAT = fem.Materials();

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

		// make sure we have a valid rigid body reference
		if ((rb < 0)||(rb>=NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

		FERigidNode* prn = new FERigidNode(&fem);

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
		else if (tag == "tol"    ) m_pim->value(tag, pLCS->m_tol);
		else if (tag == "penalty") m_pim->value(tag, pLCS->m_eps);
		else if (tag == "maxaug" ) m_pim->value(tag, pLCS->m_naugmax);
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

	int N, nf[9];

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
			else if (tag == "tri6" ) el.SetType(FE_TRI6NI );
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(m_pim->m_ntri3);
			else if (tag == "tri6" ) el.SetType(m_pim->m_ntri6);
			else if (tag == "tri7" ) el.SetType(m_pim->m_ntri7);
			else if (tag == "quad8") el.SetType(FE_QUAD8G9);
			else if (tag == "quad9") el.SetType(FE_QUAD9G9);
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

//---------------------------------------------------------------------------------
// parse a surface section for contact definitions
//
bool FEBioContactSection::BuildSurface(FESurface& s, FEFacetSet& fs, bool bnodal)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();
	int NN = m.Nodes();

	// count nr of faces
	int faces = fs.Faces();

	// allocate storage for faces
	s.create(faces);

	// read faces
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);
		FEFacetSet::FACET& fi = fs.Face(i);

		// set the element type/integration rule
		if (bnodal)
		{
			if      (fi.ntype == 4) el.SetType(FE_QUAD4NI);
			else if (fi.ntype == 3) el.SetType(FE_TRI3NI );
			else if (fi.ntype == 6) el.SetType(FE_TRI6NI );
			else return false;
		}
		else
		{
			if      (fi.ntype == 4) el.SetType(FE_QUAD4G4);
			else if (fi.ntype == 3) el.SetType(m_pim->m_ntri3);
			else if (fi.ntype == 6) el.SetType(m_pim->m_ntri6);
			else if (fi.ntype == 7) el.SetType(m_pim->m_ntri7);
			else if (fi.ntype == 8) el.SetType(FE_QUAD8G9);
			else if (fi.ntype == 9) el.SetType(FE_QUAD9G9);
			else return false;
		}

		int N = el.Nodes(); assert(N == fi.ntype);
		for (int j=0; j<N; ++j) el.m_node[j] = fi.node[j];
	}

	// copy the name
	s.SetName(fs.GetName());

	return true;
}
