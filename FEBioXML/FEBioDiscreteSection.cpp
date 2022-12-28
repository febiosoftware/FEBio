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
#include "FEBioDiscreteSection.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FECore/FEDiscreteDomain.h"
#include "FECore/FEModel.h"
#include "FECore/FEModelLoad.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
void FEBioDiscreteSection::Parse(XMLTag& tag)
{
	// make sure this tag has children
	if (tag.isleaf()) return;

	++tag;
	do
	{
		if      (tag == "spring"           ) ParseSpringSection  (tag);
		else if (tag == "rigid_axial_force") ParseRigidAxialForce(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioDiscreteSection::ParseSpringSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// determine the spring type
	const char* szt = tag.AttributeValue("type", true);

	// in 2.5 the names of the spring materials were changed so
	// we have to map the type string to the new names.
	if ((szt == 0) || (strcmp(szt, "linear") == 0)) szt = "linear spring";
	else if (strcmp(szt, "tension-only linear") == 0) szt = "tension-only linear spring";
	else if (strcmp(szt, "nonlinear") == 0) szt = "nonlinear spring";

	FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(szt, &fem));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
	spec.eclass = FE_ELEM_DISCRETE;
	spec.eshape = ET_DISCRETE;
	spec.etype  = FE_DISCRETE;
	FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(spec, &mesh, pm));
	assert(pd);
	if (pd == 0) throw FEBioImport::InvalidDomainType();
	mesh.AddDomain(pd);

	int elems = fem.GetMesh().Elements();
	int maxid = elems+1;

	// read spring discrete elements
	++tag;
	do
	{
		// read the required node tag
		if (tag == "node")
		{
			int n[2];
			tag.value(n, 2);
			n[0] -= 1;
			n[1] -= 1;
			pd->AddElement(++maxid, n);
		}
		else
		{
			// first read the domain paramters
			FEParameterList& pl = pd->GetParameterList();
			if (ReadParameter(tag, pl) == 0)
			{
				// read the actual spring material parameters
				FEParameterList& pl = pm->GetParameterList();
				if (ReadParameter(tag, pl) == 0)
				{
					throw XMLReader::InvalidTag(tag);
				}
			}
		}
		++tag;
	}
	while (!tag.isend());

	pd->CreateMaterialPointData();
}

//---------------------------------------------------------------------------------
void FEBioDiscreteSection::ParseRigidAxialForce(XMLTag& tag)
{
	// create a new rigid constraint
	FEModelLoad* paf = fecore_new<FEModelLoad>(tag.Name(), GetFEModel());

	// read the parameters
	FEParameterList& pl = paf->GetParameterList();
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add it to the model
	FEModel& fem = *GetFEModel();
	fem.AddModelLoad(paf);
}

//-----------------------------------------------------------------------------
void FEBioDiscreteSection25::Parse(XMLTag& tag)
{
	if (tag.isleaf()) return;

	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	vector<FEDiscreteMaterial*> dmat;

	int elems = fem.GetMesh().Elements();
	int maxid = elems + 1;

	++tag;
	do
	{
		if (tag == "discrete_material")
		{
			// determine the discrete material type
			const char* szt = tag.AttributeValue("type");

			// get the optional name
			const char* szname = tag.AttributeValue("name", true);

			// create the discrete material
			FEDiscreteMaterial* pm = fecore_new<FEDiscreteMaterial>(szt, &fem);
			if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

			// set the optional name
			if (szname) pm->SetName(szname);

			// add it to the model
			fem.AddMaterial(pm);

			dmat.push_back(pm);

			// read the parameter list
			ReadParameterList(tag, pm);

			// The id in the file is one-based for discrete materials, but we want to have 
			// a global material ID here.
			pm->SetID(fem.Materials());
		}
		else if (tag == "discrete")
		{
			// determine the material
			int mid;
			tag.AttributeValue("dmat", mid);
			if ((mid < 1) || (mid > (int) dmat.size())) throw XMLReader::InvalidAttributeValue(tag, "dmat");

			// create a new spring "domain"
			FE_Element_Spec spec;
			spec.eclass = FE_ELEM_DISCRETE;
			spec.eshape = ET_TRUSS2;
			spec.etype  = FE_DISCRETE;

			const char* sztype = tag.AttributeValue("type", true);
			if (sztype)
			{
				if (strcmp(sztype, "wire") == 0)
				{
					spec.eclass = FE_ELEM_WIRE;					
				}
				else if (strcmp(sztype, "spring") == 0)
				{
					spec.eclass = FE_ELEM_DISCRETE;
				}
				else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
			}

			FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(spec, &mesh, dmat[mid - 1]));
			mesh.AddDomain(pd);

			// get the discrete set
			const char* szset = tag.AttributeValue("discrete_set");
			FEDiscreteSet* pset = mesh.FindDiscreteSet(szset);
			if (pset == 0) throw XMLReader::InvalidAttributeValue(tag, "discrete_set", szset);

			// set the name of the domain
			pd->SetName(szset);

			// build the springs
			int N = pset->size();
			for (int i=0; i<N; ++i)
			{
				const FEDiscreteSet::NodePair& np = pset->Element(i);
				int n[2] = {np.n0, np.n1};
				pd->AddElement(++maxid, n);
			}

			// get the domain parameters
			FEParameterList& pl = pd->GetParameterList();
			ReadParameterList(tag, pl);

			// create an element set for this domain
			FEElementSet* elset = new FEElementSet(&fem);
			elset->Create(pd);
			elset->SetName(pd->GetName());
			mesh.AddElementSet(elset);

			// create material point data for this domain

			// initialize domain
			pd->CreateMaterialPointData();
		}
		else if (tag == "rigid_axial_force") ParseRigidAxialForce(tag);
		else if (tag == "rigid_cable") ParseRigidCable(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//---------------------------------------------------------------------------------
void FEBioDiscreteSection25::ParseRigidAxialForce(XMLTag& tag)
{
	// create a new rigid constraint
	FEModelLoad* paf = fecore_new<FEModelLoad>(tag.Name(), GetFEModel());

	// read the parameters
	FEParameterList& pl = paf->GetParameterList();
	++tag;
	do
	{
		if (ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	// add it to the model
	FEModel& fem = *GetFEModel();
	fem.AddModelLoad(paf);
}

//-----------------------------------------------------------------------------
void FEBioDiscreteSection25::ParseRigidCable(XMLTag& tag)
{
	// create a new rigid constraint
	FEModelLoad* pml = fecore_new<FEModelLoad>(tag.Name(), GetFEModel());
	if (pml == nullptr) throw XMLReader::InvalidTag(tag);

	// read all parameters and properties
	++tag;
	do
	{
		if (ReadParameter(tag, pml) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	} while (!tag.isend());

	// add it to the model
	FEModel& fem = *GetFEModel();
	fem.AddModelLoad(pml);
}
