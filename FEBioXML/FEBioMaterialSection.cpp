#include "stdafx.h"
#include "FEBioMaterialSection.h"
#include "FECore/FEModel.h"
#include "FECore/FECoreKernel.h"

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

	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// create a new material of this type
	FEMaterial* pmat = fecore_new<FEMaterial>(FEMATERIAL_ID, sztype, &fem);
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// set the material's name
	if (szname) pmat->SetName(szname);

	// set the material attributes
	for (int i=0; i<tag.m_natt; ++i)
	{
		XMLReader::XMLAtt& att = tag.m_att[i];
		if (pmat->SetAttribute(att.m_szatt, att.m_szatv) == false) { delete pmat; throw XMLReader::InvalidAttributeValue(tag, att.m_szatt); };
	}

	return pmat;
}

//-----------------------------------------------------------------------------
//! Parse the Materials section. 
//! \todo I should probably check that all tags are material tags.
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEFEBioImport::DuplicateMaterialSection();

	// reset material counter
	m_nmat = 0;

	++tag;
	do
	{
		// create a material from this tag
		FEMaterial* pmat = CreateMaterial(tag); assert(pmat);

		// add the material
		fem.AddMaterial(pmat);
		++m_nmat;

		// set the material's ID
		pmat->SetID(m_nmat);

		// parse the material parameters
		ParseMaterial(tag, pmat);

		// read next tag
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! Parse the parameters and properties of the material.
void FEBioMaterialSection::ParseMaterial(XMLTag &tag, FEMaterial* pmat)
{
	// make sure the tag is not empty
	if (tag.isleaf()) return;

	FEModel& fem = *GetFEModel();

	// get the material's parameter list
	FEParameterList& pl = pmat->GetParameterList();

	// loop over all parameters
	++tag;
	do
	{
		// see if we can find this parameter
		if (m_pim->ReadParameter(tag, pl) == false)
		{
			// if we get here the parameter was not part of the parameter list
			// however, not all parameters can be read from the parameter lists yet
			// so we have to read in the other parameters the hard way
			// TODO: this option will become obselete
			bool bfound = false;

			// additional material parameters
			if (!bfound && (tag == "fiber"   )) bfound = ParseFiberTag  (tag, pmat);
			if (!bfound && (tag == "mat_axis")) bfound = ParseMatAxisTag(tag, pmat);

			// If we get here, we use the new "material property" interface.
			if (!bfound)
			{
				// find a property with this name
				int n = pmat->FindPropertyIndex(tag.Name());
				if (n >= 0)
				{
					// create a material from this tag
					FEMaterial* pmc = CreateMaterial(tag); assert(pmc);

					// try to set the material property
					if (pmat->SetProperty(n, pmc))
					{
						bfound = true;
		
						// parse the material
						ParseMaterial(tag, pmc);
					}
				}
			}
		
			// see if we have processed the tag
			if (bfound == false) throw XMLReader::InvalidTag(tag);
		}

		// get the next tag
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
//! read the material axis
bool FEBioMaterialSection::ParseMatAxisTag(XMLTag &tag, FEMaterial *pm)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	XMLAtt& type = tag.Attribute("type");
	if (type == "local")
	{
		FELocalMap* pmap = new FELocalMap(&fem);
		pm->SetCoordinateSystemMap(pmap);

		int n[3] = {0};
		tag.value(n, 3);
		if ((n[0] == 0) && (n[1] == 0) && (n[2] == 0)) { n[0] = 1; n[1] = 2; n[2] = 4; }

		pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
	}
	else if (type == "vector")
	{
		FEVectorMap* pmap = new FEVectorMap(&fem);
		pm->SetCoordinateSystemMap(pmap);

		vec3d a(1,0,0), d(0,1,0);
		++tag;
		do
		{
			if (tag == "a") tag.value(a);
			else if (tag == "d") tag.value(d);
			else throw XMLReader::InvalidTag(tag);
			
			++tag;
		}
		while (!tag.isend());
		pmap->SetVectors(a, d);
	}
	else if (type == "user")
	{
		// material axis are read in from the ElementData section
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());

	return true;
}

//-----------------------------------------------------------------------------
//! Read the fiber tag.
bool FEBioMaterialSection::ParseFiberTag(XMLTag &tag, FEMaterial *pm)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	FECoreKernel& febio = FECoreKernel::GetInstance();

	// create a new coordinate system generator
	XMLAtt& type = tag.Attribute("type");
	if (type == "user")
	{
		fprintf(stderr, "WARNING: The ""user"" fiber type is deprecated\n");
	}
	else if (type == "local")
	{
		FELocalMap* pmap = new FELocalMap(&fem);
		pm->SetCoordinateSystemMap(pmap);

		int n[3] = {0};
		tag.value(n, 3);
		if ((n[0] == 0) && (n[1] == 0) && (n[2] == 0)) { n[0] = 1; n[1] = 2; n[2] = 4; }

		pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
	}
	else if (type == "vector")
	{
		FEVectorMap* pmap = new FEVectorMap(&fem);
		pm->SetCoordinateSystemMap(pmap);

		vec3d a;
		tag.value(a);
		pmap->SetVectors(a, vec3d(0,0,1));
	}
	else
	{
		FECoordSysMap* pmap = fecore_new<FECoordSysMap>(FECOORDSYSMAP_ID, type.cvalue(), &fem);
		if (pmap == 0) throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());

		// read the parameter list
		FEParameterList& pl = pmap->GetParameterList();
		if (pl.Parameters() > 0)
		{
			if (tag.isleaf()) m_pim->ReadParameter(tag, pl, type.cvalue());
			else
			{
				++tag;
				do
				{
					if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while(!tag.isend());
			}
		}
		pm->SetCoordinateSystemMap(pmap);
	}

	// mark the tag as read
	return true;
}
