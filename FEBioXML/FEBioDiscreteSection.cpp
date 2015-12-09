#include "stdafx.h"
#include "FEBioDiscreteSection.h"
#include "FECore/FEDiscreteMaterial.h"
#include "FECore/FEDiscreteDomain.h"
#include "FECore/FEModel.h"
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
	if (szt == 0) szt = "linear";
	FEDiscreteMaterial* pm = dynamic_cast<FEDiscreteMaterial*>(fecore_new<FEMaterial>(FEMATERIAL_ID, szt, &fem));
	if (pm == 0) throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
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

	pd->InitMaterialPointData();
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
