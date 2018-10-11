#include "stdafx.h"
#include "FEBioMaterialSection.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! This function creates a material by checking the type attribute against
//! registered materials. Also, if the tag defines attributes (other than
//! type and name), the material is offered a chance to process the attributes.
FEMaterial* FEBioMaterialSection::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);
	
	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// create a new material of this type
	FEMaterial* pmat = GetBuilder()->CreateMaterial(sztype);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	return pmat;
}

//-----------------------------------------------------------------------------
// This helper function checks if all the required material properties are defined
void ValidateMaterial(FEMaterial* pmat)
{
	int NP = pmat->PropertyClasses();
	for (int i=0; i<NP; ++i)
	{
		FEProperty* pi = pmat->PropertyClass(i);
		if (pi->m_brequired && (pi->size()==0))
		{
			throw FEBioImport::MissingMaterialProperty(pmat->GetName(), pi->GetName());
		}
	}
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		if (tag == "material")
		{
			// create a material from this tag
			FEMaterial* pmat = CreateMaterial(tag); assert(pmat);

			// check that the ID attribute is defined and that it 
			// equals the number of materials + 1.
			int nid = -1;
			tag.AttributeValue("id", nid);
			int nmat = fem.Materials();
			if (nid != nmat+1) throw XMLReader::InvalidAttributeValue(tag, "id");

			// add the material
			fem.AddMaterial(pmat);
			++m_nmat;

			// set the material's ID
			pmat->SetID(m_nmat);

			// parse the material parameters
			ReadParameterList(tag, pmat);

			// validate the material
			ValidateMaterial(pmat);
		}
		else throw XMLReader::InvalidTag(tag);

		// read next tag
		++tag;
	}
	while (!tag.isend());
}
