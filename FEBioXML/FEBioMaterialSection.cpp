#include "stdafx.h"
#include "FEBioMaterialSection.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
FEMaterial* FEBioMaterialSection::CreateMaterial(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	FEBioKernel& febio = FEBioKernel::GetInstance();

	// get the material type
	const char* sztype = tag.AttributeValue("type", true);
	
	// in some case, a type is not defined (e.g. for solutes)
	// in that case, we use the tag name as the type
	if (sztype == 0) sztype = tag.Name();

	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// create a new material of this type
	FEMaterial* pmat = febio.Create<FEMaterial>(sztype, &fem);
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
void FEBioMaterialSection::ParseMaterial(XMLTag &tag, FEMaterial* pmat)
{
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

			// additional elastic material parameters
			if (!bfound && dynamic_cast<FEElasticMaterial*>(pmat)) bfound = ParseElasticMaterial(tag, dynamic_cast<FEElasticMaterial*>(pmat));

			// additional transversely isotropic material parameters
			if (!bfound && dynamic_cast<FETransverselyIsotropic*>(pmat)) bfound = ParseTransIsoMaterial(tag, dynamic_cast<FETransverselyIsotropic*>(pmat));

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
// Parse FEElasticMaterial 
//
bool FEBioMaterialSection::ParseElasticMaterial(XMLTag &tag, FEElasticMaterial *pm)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	// read the material axis
	if (tag == "mat_axis")
	{
		XMLAtt& type = tag.Attribute("type");
		if (type == "local")
		{
			FELocalMap* pmap = new FELocalMap(&fem);
			pm->m_pmap = pmap;

			int n[3] = {0};
			tag.value(n, 3);
			if ((n[0] == 0) && (n[1] == 0) && (n[2] == 0)) { n[0] = 1; n[1] = 2; n[2] = 4; }

			pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
		}
		else if (type == "vector")
		{
			FEVectorMap* pmap = new FEVectorMap(&fem);
			pm->m_pmap = pmap;

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
	else if (tag == "fiber")
	{
		FEBioKernel& febio = FEBioKernel::GetInstance();

		// create a new coordinate system generator
		XMLAtt& type = tag.Attribute("type");
		if (type == "user")
		{
			fprintf(stderr, "WARNING: The ""user"" fiber type is deprecated\n");
		}
		else
		{
			FECoordSysMap* pmap = febio.Create<FECoordSysMap>(type.cvalue(), &fem);
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
			pm->m_pmap = pmap;
		}

		// mark the tag as read
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Parse FETransverselyIsotropic
//
bool FEBioMaterialSection::ParseTransIsoMaterial(XMLTag &tag, FETransverselyIsotropic *pm)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	if (tag == "active_contraction")
	{
		const char* szlc = tag.AttributeValue("lc");
		FEParameterList& pl = pm->m_fib.GetParameterList();
		FEParam& p = *pl.Find("ascl");
		p.m_nlc = atoi(szlc)-1;
		p.value<double>() = 1.0;

		++tag;
		do
		{
			if (m_pim->ReadParameter(tag, pl) == false) throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());

		// mark tag as read
		return true;
	}
	return false;
}
