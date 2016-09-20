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
	int nversion = m_pim->Version();
	if (nversion < 0x0205) ParseDiscreteSection(tag);
	else ParseDiscreteSection25(tag);
}

//-----------------------------------------------------------------------------
void FEBioDiscreteSection::ParseDiscreteSection(XMLTag& tag)
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

	FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, szt, &fem));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
	spec.eclass = FE_ELEM_TRUSS;
	spec.eshape = ET_TRUSS2;
	spec.etype  = FE_DISCRETE;
	FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(spec, &mesh, pm));
	mesh.AddDomain(pd);

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
			pd->AddElement(++m_pim->m_maxid, n);
		}
		else
		{
			// first read the domain paramters
			FEParameterList& pl = pd->GetParameterList();
			if (m_pim->ReadParameter(tag, pl) == 0)
			{
				// read the actual spring material parameters
				FEParameterList& pl = pm->GetParameterList();
				if (m_pim->ReadParameter(tag, pl) == 0)
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
	FEModelLoad* paf = fecore_new<FEModelLoad>(FEBC_ID, tag.Name(), GetFEModel());

	// read the parameters
	FEParameterList& pl = paf->GetParameterList();
	++tag;
	do
	{
		if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add it to the model
	FEModel& fem = *GetFEModel();
	fem.AddModelLoad(paf);
}

//-----------------------------------------------------------------------------
void FEBioDiscreteSection::ParseDiscreteSection25(XMLTag& tag)
{
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	vector<FEDiscreteMaterial*> dmat;

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
			FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, szt, &fem));
			if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

			// set the optional name
			if (szname) pm->SetName(szname);

			// add it to the model
			fem.AddMaterial(pm);
			pm->SetID(fem.Materials());

			dmat.push_back(pm);

			// read the parameter list
			FEParameterList& pl = pm->GetParameterList();
			++tag;
			do
			{
				if (m_pim->ReadParameter(tag, pl) == 0) throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
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

			// build the springs
			int N = pset->size();
			for (int i=0; i<N; ++i)
			{
				const FEDiscreteSet::NodePair& np = pset->Element(i);
				int n[2] = {np.n0, np.n1};
				pd->AddElement(++m_pim->m_maxid, n);
			}

			// get the domain parameters
			FEParameterList& pl = pd->GetParameterList();
			m_pim->ReadParameterList(tag, pl);

			// initialize domain
			pd->CreateMaterialPointData();
		}
		else if (tag == "rigid_axial_force") ParseRigidAxialForce(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
