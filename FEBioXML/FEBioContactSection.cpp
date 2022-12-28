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
#include "FEBioContactSection.h"
#include <FECore/FEAugLagLinearConstraint.h>
#include <FECore/FENodalBC.h>
#include <FECore/FECoreKernel.h>
#include <FECore/FEModel.h>

//-----------------------------------------------------------------------------
FEBioContactSection::MissingPrimarySurface::MissingPrimarySurface()
{
	SetErrorString("Missing contact primary surface");
}

//-----------------------------------------------------------------------------
FEBioContactSection::MissingSecondarySurface::MissingSecondarySurface()
{
	SetErrorString("Missing contact secondary surface");
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
					else
					{
						// check for some older obsolete names
						if      (strcmp(sztype, "facet-to-facet sliding") == 0) pci = fecore_new_class<FESurfacePairConstraint>("FEFacet2FacetSliding", &fem);
						else if (strcmp(sztype, "sliding_with_gaps"     ) == 0) pci = fecore_new_class<FESurfacePairConstraint>("FESlidingInterface", &fem);

						if (pci)
						{
							GetBuilder()->AddContactInterface(pci);
							ParseContactInterface(tag, pci);
						}
						else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
					}
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

				if (pci == nullptr)
				{
					// check for some older obsolete names
					if (strcmp(sztype, "facet-to-facet sliding") == 0) pci = fecore_new_class<FESurfacePairConstraint>("FEFacet2FacetSliding", &fem);
					else if (strcmp(sztype, "sliding_with_gaps") == 0) pci = fecore_new_class<FESurfacePairConstraint>("FESlidingInterface", &fem);
				}

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

				FESurface& s = *(ntype == 1? pci->GetSecondarySurface() : pci->GetPrimarySurface());
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

	// Make sure we have a primary and secondary interface
	FESurface* pss = pci->GetPrimarySurface (); if ((pss == 0) || (pss->Elements()==0)) throw MissingPrimarySurface ();
	FESurface* pms = pci->GetSecondarySurface(); if ((pms == 0) || (pms->Elements()==0)) throw MissingSecondarySurface();
}

//-----------------------------------------------------------------------------
void FEBioContactSection25::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the surface pair
	const char* szpair = tag.AttributeValue("surface_pair");
	FESurfacePair* surfacePair =m.FindSurfacePair(szpair);
	if (surfacePair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// build the surfaces
	if (GetBuilder()->BuildSurface(*pci->GetSecondarySurface(), *surfacePair->GetSecondarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
	if (GetBuilder()->BuildSurface(*pci->GetPrimarySurface(), *surfacePair->GetPrimarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// get the parameter list
	ReadParameterList(tag, pci);

	// Make sure we have both surfaces
	FESurface* pss = pci->GetPrimarySurface (); if ((pss == 0) || (pss->Elements()==0)) throw MissingPrimarySurface ();
	m.AddSurface(pss);
	FESurface* pms = pci->GetSecondarySurface(); if ((pms == 0) || (pms->Elements()==0)) throw MissingSecondarySurface();
	m.AddSurface(pms);
}

//-----------------------------------------------------------------------------
// --- R I G I D   W A L L   I N T E R F A C E ---
void FEBioContactSection2::ParseRigidWall(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	FESurfacePairConstraint* ps = fecore_new<FESurfacePairConstraint>("rigid_wall", GetFEModel());
	fem.AddSurfacePairConstraint(ps);

	++tag;
	do
	{
		if (ReadParameter(tag, ps) == false)
		{
			if (tag == "surface")
			{
				FESurface& s = *ps->GetPrimarySurface();

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if      (strcmp(szfmt, "face nodes"  ) == 0) nfmt = 0;
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

	FESurfaceConstraint* ps = fecore_new<FESurfaceConstraint>("rigid_wall", GetFEModel());
	fem.AddNonlinearConstraint(ps);

	// get and build the surface
	const char* sz = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(sz);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	if (GetBuilder()->BuildSurface(*ps->GetSurface(), *pface, true) == false) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);

	mesh.AddSurface(ps->GetSurface());

	ReadParameterList(tag, ps);
}

//-----------------------------------------------------------------------------
void FEBioContactSection25::ParseRigidSliding(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FESurfacePairConstraint* ps = fecore_new<FESurfacePairConstraint>("rigid sliding", GetFEModel());
	fem.AddSurfacePairConstraint(ps);

	// get and build the surface
	const char* sz = tag.AttributeValue("surface");
	FEFacetSet* pface = mesh.FindFacetSet(sz);
	if (pface == 0) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	if (GetBuilder()->BuildSurface(*ps->GetPrimarySurface(), *pface, false) == false) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);
	mesh.AddSurface(ps->GetPrimarySurface());

	ReadParameterList(tag, ps);
}

//-----------------------------------------------------------------------------
// --- R I G I D   B O D Y   I N T E R F A C E ---
void FEBioContactSection2::ParseRigidInterface(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEModelBuilder* feb = GetBuilder();

	int NMAT = fem.Materials();

	// count how many rigid nodes there are
	int nrn= 0;
	XMLTag t(tag); ++t;
	while (!t.isend()) { nrn++; ++t; }

	++tag;
	int id, rb, rbp = -1;
	FENodalBC* prn = 0;
	FENodeSet* ns = nullptr;
	for (int i=0; i<nrn; ++i)
	{
		id = atoi(tag.AttributeValue("id"))-1;
		rb = atoi(tag.AttributeValue("rb"));

		// make sure we have a valid rigid body reference
		if ((rb <= 0)||(rb>NMAT)) throw XMLReader::InvalidAttributeValue(tag, "rb", tag.AttributeValue("rb"));

		if ((prn == 0) || (rb != rbp))
		{
			prn = fecore_new_class<FENodalBC>("FERigidNodeSet", &fem);
			ns = new FENodeSet(&fem);
			prn->SetNodeSet(ns);

			// the default shell bc depends on the shell formulation
			// hinged shell = 0
			// clamped shell = 1
			prn->SetParameter("clamp_shells", feb->m_default_shell == OLD_SHELL ? 0 : 1);

			prn->SetParameter("rb", rb);

			GetBuilder()->AddBC(prn);
			rbp = rb;
		}
		if (ns) ns->Add(id);

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
			FEAugLagLinearConstraint* pLC = fecore_alloc(FEAugLagLinearConstraint, &fem);

			++tag;
			do
			{
				int node, bc;
				double val;
				if (tag == "node")
				{
					tag.value(val);
					
					const char* szid = tag.AttributeValue("id");
					node = atoi(szid);

					const char* szbc = tag.AttributeValue("bc");
                    int ndof = dofs.GetDOF(szbc);
                    if (ndof >= 0) bc = ndof;
					else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

					pLC->AddDOF(node, bc, val);
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
				int nn = pe->GetFace(nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
				el.m_elem[0] = pe;
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}

	s.InitSurface();
	s.CreateMaterialPointData();

	return true;
}
//-----------------------------------------------------------------------------
//! Parse the Contact section (new in version 2.0)
void FEBioContactSection4::Parse(XMLTag& tag)
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

			// Try to create a contact interface using the FEBio kernel. 
			FESurfacePairConstraint* pci = fecore_new<FESurfacePairConstraint>(sztype, &fem);
			if (pci == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			GetBuilder()->AddContactInterface(pci);
			ParseContactInterface(tag, pci);
		}

		++tag;
	} while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioContactSection4::ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci)
{
	FEModel& fem = *GetFEModel();
	FEMesh& m = fem.GetMesh();

	// get the surface pair
	const char* szpair = tag.AttributeValue("surface_pair");
	FESurfacePair* surfacePair = m.FindSurfacePair(szpair);
	if (surfacePair == 0) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// build the surfaces
	if (GetBuilder()->BuildSurface(*pci->GetSecondarySurface(), *surfacePair->GetSecondarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);
	if (GetBuilder()->BuildSurface(*pci->GetPrimarySurface(), *surfacePair->GetPrimarySurface(), pci->UseNodalIntegration()) == false) throw XMLReader::InvalidAttributeValue(tag, "surface_pair", szpair);

	// get the parameter list
	ReadParameterList(tag, pci);

	// Make sure we have both surfaces
	FESurface* pss = pci->GetPrimarySurface(); if ((pss == 0) || (pss->Elements() == 0)) throw MissingPrimarySurface();
	m.AddSurface(pss);
	FESurface* pms = pci->GetSecondarySurface(); if ((pms == 0) || (pms->Elements() == 0)) throw MissingSecondarySurface();
	m.AddSurface(pms);
}
