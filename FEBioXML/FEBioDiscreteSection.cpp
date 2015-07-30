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

	// create a new spring "domain"
	FECoreKernel& febio = FECoreKernel::GetInstance();
	FE_Element_Spec spec;
	spec.eshape = ET_TRUSS2;
	spec.etype  = FE_DISCRETE;
	int ndomtype = febio.GetDomainType(spec, pm);
	FEDiscreteDomain* pd = dynamic_cast<FEDiscreteDomain*>(febio.CreateDomain(ndomtype, &mesh, pm));
	mesh.AddDomain(pd);

	pd->create(1);
	FEDiscreteElement& de = pd->Element(0);
	de.SetType(FE_DISCRETE);
	de.m_nID = ++m_pim->m_maxid;
	
	// add a new material for each spring
	fem.AddMaterial(pm);
	pm->SetID(fem.Materials());
	de.SetMatID(fem.Materials()-1);

	// read spring discrete elements
	++tag;
	do
	{
		// read the required node tag
		if (tag == "node")
		{
			int n[2];
			tag.value(n, 2);
			de.m_node[0] = n[0]-1;
			de.m_node[1] = n[1]-1;
		}
		else
		{
			// read the actual spring material parameters
			FEParameterList& pl = pm->GetParameterList();
			if (m_pim->ReadParameter(tag, pl) == 0)
			{
				throw XMLReader::InvalidTag(tag);
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
	FERigidAxialForce* paf = new FERigidAxialForce(GetFEModel());

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
	fem.m_RAF.push_back(paf);
}
