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

	// TODO: Solute materials need to have their ID set
	if (dynamic_cast<FESolute*>(pmat))
	{
		FESolute* ps = dynamic_cast<FESolute*>(pmat);
		const char* szid = tag.AttributeValue("sol");
		int nid = atoi(szid) - 1;
		if ((nid < 0) || (nid >= MAX_CDOFS)) throw XMLReader::InvalidAttributeValue(tag, "sol", szid);
		ps->SetSoluteID(nid);
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

			// chemical reaction material parameters
			if (!bfound && dynamic_cast<FEChemicalReaction*>(pmat)) bfound = ParseReactionMaterial(tag, dynamic_cast<FEChemicalReaction*>(pmat));
			
			// multiphasic material parameters
			if (!bfound && dynamic_cast<FEMultiphasic*>(pmat)) bfound = ParseMultiphasicMaterial(tag, dynamic_cast<FEMultiphasic*>(pmat));

			// multigeneration materials
			if (!bfound && dynamic_cast<FEElasticMultigeneration*>(pmat)) bfound = ParseElasticMultigeneration(tag, dynamic_cast<FEElasticMultigeneration*>(pmat));

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

	// read material fibers
	if (tag == "fiber")
	{
		FEBioKernel& febio = FEBioKernel::GetInstance();

		// create a new coordinate system generator
		const char* sztype = tag.AttributeValue("type");
		if (strcmp(sztype, "user") == 0)
		{
			fprintf(stderr, "WARNING: The ""user"" fiber type is deprecated\n");
		}
		else
		{
			FECoordSysMap* pmap = febio.Create<FECoordSysMap>(sztype, &fem);
			if (pmap == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			// read the parameter list
			FEParameterList& pl = pmap->GetParameterList();
			if (pl.Parameters() > 0)
			{
				if (tag.isleaf()) m_pim->ReadParameter(tag, pl, sztype);
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
	else if (tag == "active_contraction")
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

//-----------------------------------------------------------------------------
// Parse ParseElasticMultigeneration material 
//
bool FEBioMaterialSection::ParseElasticMultigeneration(XMLTag &tag, FEElasticMultigeneration *pm)
{
	// read the solid material
	if (tag == "solid")
	{
		// create a new material of this type
		FEMaterial* pmat = CreateMaterial(tag);
		
		// make sure the base material is a valid material (i.e. an elastic material)
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || pme->IsRigid()) throw XMLReader::InvalidTag(tag);
		
		// get the growth generation and store it in the material ID
		int id;
		tag.AttributeValue("gen", id);

		// make the ID zero-based
		id--;
		if (id < 0) throw XMLReader::InvalidAttributeValue(tag, "gen");

		// Set the ID
		pme->SetID(id);
		
		// add the material as a new generation
		pm->AddMaterial(pme);
		
		// parse the solid
		ParseMaterial(tag, pmat);
		
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Parse FEChemicalReaction material 
//
bool FEBioMaterialSection::ParseReactionMaterial(XMLTag &tag, FEChemicalReaction *pm)
{
	// read the stoichiometric coefficients and reaction rates
	if (tag == "vR")
	{
		int id, vR;
		if (tag.AttributeValue("sbm", id, true))
		{
			if (id < 1) throw XMLReader::InvalidAttributeValue(tag, "sbm");
			tag.value(vR);
			pm->SetSolidReactantsCoefficients(id-1, vR);
		}
		if (tag.AttributeValue("sol", id, true))
		{
			if ((id < 1) || (id > MAX_CDOFS)) throw XMLReader::InvalidAttributeValue(tag, "sol" );
			tag.value(vR);
			pm->SetSoluteReactantsCoefficients(id-1, vR);
		}
		return true;
	}
	else if (tag == "vP")
	{
		int id, vP;
		if (tag.AttributeValue("sbm", id, true))
		{
			if (id < 1) throw XMLReader::InvalidAttributeValue(tag, "sbm");
			tag.value(vP);
			pm->SetSolidProductsCoefficients(id-1, vP);
		}
		if (tag.AttributeValue("sol", id, true))
		{
			if ((id < 1) || (id > MAX_CDOFS)) throw XMLReader::InvalidAttributeValue(tag, "sol");
			tag.value(vP);
			pm->SetSoluteProductsCoefficients(id-1, vP);
		}
		
		return true;
	}
	else if (tag == "forward_rate")
	{
		// create a new material of this type
		FEMaterial* pmat = CreateMaterial(tag);
		
		// make sure the base material is a valid material
		FEReactionRate* pme = dynamic_cast<FEReactionRate*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid forward reaction rate in multiphasic material %s\n", pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the reaction rate pointer
		pm->SetForwardReactionRate(pme);
		
		// parse the material
		ParseMaterial(tag, pme);
		
		return true;
	}
	else if (tag == "reverse_rate")
	{
		// create a new material of this type
		FEMaterial* pmat = CreateMaterial(tag);
		
		// make sure the base material is a valid material
		FEReactionRate* pme = dynamic_cast<FEReactionRate*>(pmat);
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid reverse reaction rate in multiphasic material %s\n", pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the reaction rate pointer
		pm->SetReverseReactionRate(pme);
		
		// parse the material
		ParseMaterial(tag, pme);
		
		return true;
	}
	else if (tag == "Vbar")
    {
        pm->m_Vovr = true;
        tag.value(pm->m_Vbar);
		return true;
    }
	else throw XMLReader::InvalidTag(tag);
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse FETriphasic material 
//
bool FEBioMaterialSection::ParseMultiphasicMaterial(XMLTag &tag, FEMultiphasic *pm)
{
	// solutes are handled slightly differently
	if (tag == "solute")
	{
		// get the (optional) material name
		const char* szname = tag.AttributeValue("name", true);

		// get the solute ID
		const char* szid = tag.AttributeValue("sol");
		int nid = atoi(szid) - 1;

		// create a new solute material
		FESolute* psol = new FESolute;
		psol->SetSoluteID(nid);

		// add the solute
		pm->AddSolute(psol);

		// set the new material's name (if defined)
		if (szname) psol->SetName(szname);
		
		// parse the new material
		ParseMaterial(tag, psol);
	}
	else if (tag == "solid_bound")
	{
		// get the material type
		const char* sztype = "solid_bound";
		
		// get the material name
		const char* szname = tag.AttributeValue("name", true);
		
		// get the solid-bound molecule id
		const char* szid = tag.AttributeValue("sbm");
		int id = atoi(szid) - 1;
		if (id < 0) throw XMLReader::InvalidAttributeValue(tag, "sbm", szid);
		
		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. a solid-bound molecule material)
		FESolidBoundMolecule* pme = dynamic_cast<FESolidBoundMolecule*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid solid-bound molecule %s in multiphasic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid-bound molecule pointer
		pm->AddSolidBoundMolecule(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// set the solute ID
		pme->SetSBMID(id);
		
		// parse the material
		ParseMaterial(tag, pme);
	}
	else if (tag == "reaction")
	{
		// create a new material of this type
		FEMaterial* pmat = CreateMaterial(tag);
		
		// make sure the base material is a valid material (i.e. a chemical reaction material)
		FEChemicalReaction* pme = dynamic_cast<FEChemicalReaction*>(pmat);
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid chemical reaction in multiphasic material %s\n", pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the chemical reaction pointer
		pm->AddChemicalReaction(pme);
		
		// parse the material
		ParseMaterial(tag, pme);
	}
	else
	{
		// see if we can find a material property with this name
		int nc = pm->FindComponent(tag.Name());
		if (nc == -1) throw XMLReader::InvalidTag(tag);

		// create a new material of this type
		FEMaterial* pmat = CreateMaterial(tag);

		// assign the new material to the corresponding material property
		if (pm->SetComponent(nc, pmat) == false)
		{
			clog.printbox("INPUT ERROR: Invalid %s definition in material %s\n", tag.Name(), pm->GetName());
			throw XMLReader::Error();
		}

		// parse the new material
		ParseMaterial(tag, pmat);
	}
	
	return true;
}
