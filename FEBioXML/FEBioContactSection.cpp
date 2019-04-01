/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEBioContactSection.h"
#include "FEBioMech/FERigidWallInterface.h"
#include "FEBioMech/FEAugLagLinearConstraint.h"
#include <FEBioMech/FERigidSlidingContact.h>
#include "FECore/FECoreKernel.h"
#include <FECore/FEModel.h>
#include <FEBioMech/FERigidSystem.h>
#include <FEBioMech/RigidBC.h>
#include <FEBioMech/FEMechModel.h>

//-----------------------------------------------------------------------------
FEBioContactSection::MissingSlaveSurface::MissingSlaveSurface()
{
	SetErrorString("Missing contact slave surface");
}

//-----------------------------------------------------------------------------
FEBioContactSection::MissingMasterSurface::MissingMasterSurface()
{
	SetErrorString("Missing contact master surface");
}

//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection2::Parse(XMLTag& tag)
{
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
			if      (strcmp(sztype, "rigid_wall"       ) == 0) ParseRigidWall       (tag);
			else if (strcmp(sztype, "rigid"            ) == 0) ParseRigidInterface  (tag);
			else if (strcmp(sztype, "linear constraint") == 0) ParseLinearConstraint(tag);
			else 
			{
				// If we get here, we try to create a contact interface
				// using the FEBio kernel. 
				FESurfacePairConstraint* pci = fecore_new<FESurfacePairConstraint>(sztype, &fem);
				if (pci)
				{
					GetBuilder()->AddContactInterface(pci);
					ParseContactInterface(tag, pci);
				}
				else
				{
					FEModel& fem = *GetFEModel();

					// Some constraints were initially defined in the Contact section, although
					// now it is preferred that they are defined in the Constraints section. For backward
					// compatibility we still allow constraints to be defined in this section. 
					FENLConstraint* pc = fecore_new<FENLConstraint>(sztype, &fem);
					if (pc)
					{
						ReadParameterList(tag, pc);
						GetBuilder()->AddNonlinearConstraint(pc);
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
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection25::Parse(XMLTag& tag)
{
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
			if      (strcmp(sztype, "rigid_wall"       ) == 0) ParseRigidWall       (tag);
			else if (strcmp(sztype, "rigid sliding"    ) == 0) ParseRigidSliding    (tag);
			else if (strcmp(sztype, "linear constraint") == 0) ParseLinearConstraint(tag);
			else
			{
				// If we get here, we try to create a contact interface
				// using the FEBio kernel. 
				FESurfacePairConstraint* pci = fecore_new<FESurfacePairConstraint>(sztype, &fem);
				if (pci)
				{
					GetBuilder()->AddContactInterface(pci);
					ParseContactInterface(tag, pci);
				}
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioContactSection2::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();

	// read the parameters
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false)
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
					FEFacetSet* ps = m.FindFacetSet(szset);

					// create a surface from the facet set
					if (ps)
					{
						if (GetBuilder()->BuildSurface(s, *ps, pci->UseNodalIntegration()) == false) throw XMLReader::InvalidTag(tag);
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

	// Make sure we have a master and a slave interface
	FESurface* pss = pci->GetSlaveSurface (); if ((pss == 0) || (pss->Elements()==0)) throw MissingSlaveSurface ();
	FESurface* pms = pci->GetMasterSurface(); if ((pms == 0) || (pms->Elements()==0)) throw MissingMasterSurface();
}

//-----------------------------------------------------------------------------
void FEBioContactSection25::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the surface pair
	const char* szpair = tag.AttributeValue("surface_pair");
	FESurfacePair* surfacePair = m.FindSurfacePair(szpair);
	if (surfacePair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// build the surfaces
	if (GetBuilder()->BuildSurface(*pci->GetMasterSurface(), *surfacePair->GetMasterSurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
	if (GetBuilder()->BuildSurface(*pci->GetSlaveSurface(), *surfacePair->GetSlaveSurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// get the parameter list
	FEParameterList& pl = pci->GetParameterList();
	ReadParameterList(tag, pl);

	// Make sure we have a master and a slave interface
	FESurface* pss = pci->GetSlaveSurface (); if ((pss == 0) || (pss->Elements()==0)) throw MissingSlaveSurface ();
	m.AddSurface(pss);
	FESurface* pms = pci->GetMasterSurface(); if ((pms == 0) || (pms->Elements()==0)) throw MissingMasterSurface();
	m.AddSurface(pms);
}

//-----------------------------------------------------------------------------
// --- R I G I D   W A L L   I N T E R F A C E ---
void FEBioContactSection2::ParseRigidWall(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	FERigidWallInterface* ps = dynamic_cast<FERigidWallInterface*>(fecore_new<FESurfacePairConstraint>("rigid_wall", GetFEModel()));
	fem.AddSurfacePairConstraint(ps);

	++tag;
	do
	{
		if (ReadParameter(tag, ps) == false)
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
// --- R I G I D   W A L L   I N T E R F A C E ---
void FEBioContactSection25::ParseRigidWall(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FERigidWallInterface* ps = dynamic_cast<FERigidWallInterface*>(fecore_new<FESurfacePairConstraint>("rigid_wall", GetFEModel()));
	fem.AddSurfacePairConstraint(ps);

	// get and build the surface
	const char* sz = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(sz);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	if (GetBuilder()->BuildSurface(ps->m_ss, *pface, ps->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);

	ReadParameterList(tag, ps);
}

//-----------------------------------------------------------------------------
void FEBioContactSection25::ParseRigidSliding(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FERigidSlidingContact* ps = fecore_new<FERigidSlidingContact>("rigid sliding", GetFEModel());
	fem.AddSurfacePairConstraint(ps);

	// get and build the surface
	const char* sz = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(sz);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	if (GetBuilder()->BuildSurface(*ps->GetSlaveSurface(), *pface, false) == false) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	mesh.AddSurface(ps->GetSlaveSurface());

	ReadParameterList(tag, ps);
}

//-----------------------------------------------------------------------------
// --- R I G I D   B O D Y   I N T E R F A C E ---
void FEBioContactSection2::ParseRigidInterface(XMLTag& tag)
{
	FEMechModel& fem = static_cast<FEMechModel&>(*GetFEModel());
	FERigidSystem& rigid = *fem.GetRigidSystem();
	FEModelBuilder* feb = GetBuilder();

	int NMAT = fem.Materials();

	// count how many rigid nodes there are
	int nrn= 0;
	XMLTag t(tag); ++t;
	while (!t.isend()) { nrn++; ++t; }

	++tag;
	int id, rb, rbp = -1;
	FERigidNodeSet* prn = 0;
	for (int i=0; i<nrn; ++i)
	{
		id = atoi(tag.AttributeValue("id"))-1;
		rb = atoi(tag.AttributeValue("rb"))-1;

		// make sure we have a valid rigid body reference
		if ((rb < 0)||(rb>=NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

		if ((prn == 0) || (rb != rbp))
		{
			prn = new FERigidNodeSet(&fem);

			// the default shell bc depends on the shell formulation
			prn->SetShellBC(feb->m_default_shell == OLD_SHELL ? FERigidNodeSet::HINGED_SHELL : FERigidNodeSet::CLAMPED_SHELL);

			prn->SetRigidID(rb);

			GetBuilder()->AddRigidNodeSet(prn);
			rbp = rb;
		}
		prn->AddNode(id);

		++tag;
	}
}

//-----------------------------------------------------------------------------
// --- L I N E A R   C O N S T R A I N T ---
void FEBioContactSection::ParseLinearConstraint(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
    DOFS& dofs = fem.GetDOFS();
	FEMesh& m = fem.GetMesh();

	// make sure there is a constraint defined
	if (tag.isleaf()) return;

	// create a new linear constraint manager
	FELinearConstraintSet* pLCS = dynamic_cast<FELinearConstraintSet*>(fecore_new<FENLConstraint>("linear constraint", GetFEModel()));
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
					dof.node = ReadNodeID(tag);

					const char* szbc = tag.AttributeValue("bc");
                    int ndof = dofs.GetDOF(szbc);
                    if (ndof >= 0) dof.bc = ndof;
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
		else if (ReadParameter(tag, pLCS) == false)
		{
			throw XMLReader::InvalidTag(tag);
		}
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
	s.Create(faces);

	FEModelBuilder* feb = GetBuilder();

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
            else if (tag == "quad8" ) el.SetType(FE_QUAD8NI);
            else if (tag == "quad9" ) el.SetType(FE_QUAD9NI);
			else throw XMLReader::InvalidTag(tag);
		}
		else
		{
			if      (tag == "quad4") el.SetType(FE_QUAD4G4);
			else if (tag == "tri3" ) el.SetType(feb->m_ntri3);
			else if (tag == "tri6" ) el.SetType(feb->m_ntri6);
			else if (tag == "tri7" ) el.SetType(feb->m_ntri7);
			else if (tag == "tri10") el.SetType(feb->m_ntri10);
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
				el.m_elem[0] = pe;
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}
