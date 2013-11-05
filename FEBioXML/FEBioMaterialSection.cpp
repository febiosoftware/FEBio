#include "stdafx.h"
#include "FEBioMaterialSection.h"
#include "FECore/FEModel.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEFEBioImport::DuplicateMaterialSection();

	m_nmat = 0;

	FEBioKernel& febio = FEBioKernel::GetInstance();

	++tag;
	do
	{
		// get the material type
		const char* sztype = tag.AttributeValue("type");

		// get the material name
		const char* szname = tag.AttributeValue("name", true);

		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, &fem);
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// add the material
		fem.AddMaterial(pmat);
		++m_nmat;

		// set the material's name
		if (szname) pmat->SetName(szname);

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

			// additional remodeling solid
			if (!bfound && dynamic_cast<FERemodelingElasticMaterial*>(pmat)) bfound = ParseRemodelingSolid(tag, dynamic_cast<FERemodelingElasticMaterial*>(pmat));

			// additional elastic material parameters
			if (!bfound && dynamic_cast<FEElasticMaterial*>(pmat)) bfound = ParseElasticMaterial(tag, dynamic_cast<FEElasticMaterial*>(pmat));

			// additional transversely isotropic material parameters
			if (!bfound && dynamic_cast<FETransverselyIsotropic*>(pmat)) bfound = ParseTransIsoMaterial(tag, dynamic_cast<FETransverselyIsotropic*>(pmat));

			// read rigid body data
			if (!bfound && dynamic_cast<FERigidMaterial*>(pmat)) bfound = ParseRigidMaterial(tag, dynamic_cast<FERigidMaterial*>(pmat));

			// elastic mixtures
			if (!bfound && dynamic_cast<FEElasticMixture*>(pmat)) bfound = ParseElasticMixture(tag, dynamic_cast<FEElasticMixture*>(pmat));
			
			// uncoupled elastic mixtures
			if (!bfound && dynamic_cast<FEUncoupledElasticMixture*>(pmat)) bfound = ParseUncoupledElasticMixture(tag, dynamic_cast<FEUncoupledElasticMixture*>(pmat));
			
			// biphasic material parameters
			if (!bfound && dynamic_cast<FEBiphasic*>(pmat)) bfound = ParseBiphasicMaterial(tag, dynamic_cast<FEBiphasic*>(pmat));
			
			// biphasic-solute material parameters
			if (!bfound && dynamic_cast<FEBiphasicSolute*>(pmat)) bfound = ParseBiphasicSoluteMaterial(tag, dynamic_cast<FEBiphasicSolute*>(pmat));

			// solute material parameters
			if (!bfound && dynamic_cast<FESolute*>(pmat)) bfound = ParseSoluteMaterial(tag, dynamic_cast<FESolute*>(pmat));

			// chemical reaction material parameters
			if (!bfound && dynamic_cast<FEChemicalReaction*>(pmat)) bfound = ParseReactionMaterial(tag, dynamic_cast<FEChemicalReaction*>(pmat));
			
			// triphasic material parameters
			if (!bfound && dynamic_cast<FETriphasic*>(pmat)) bfound = ParseTriphasicMaterial(tag, dynamic_cast<FETriphasic*>(pmat));

			// multiphasic material parameters
			if (!bfound && dynamic_cast<FEMultiphasic*>(pmat)) bfound = ParseMultiphasicMaterial(tag, dynamic_cast<FEMultiphasic*>(pmat));

			// viscoelastic materials
			if (!bfound && dynamic_cast<FEViscoElasticMaterial*>(pmat)) bfound = ParseViscoElasticMaterial(tag, dynamic_cast<FEViscoElasticMaterial*>(pmat));

			// uncoupled viscoelastic materials
			if (!bfound && dynamic_cast<FEUncoupledViscoElasticMaterial*>(pmat)) bfound = ParseUncoupledViscoElasticMaterial(tag, dynamic_cast<FEUncoupledViscoElasticMaterial*>(pmat));
			
			// multigeneration materials
			if (!bfound && dynamic_cast<FEElasticMultigeneration*>(pmat)) bfound = ParseElasticMultigeneration(tag, dynamic_cast<FEElasticMultigeneration*>(pmat));

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
			if      (tag == "Tmax" ) tag.value(pm->m_fib.m_Tmax );
			else if (tag == "ca0"  ) tag.value(pm->m_fib.m_ca0  );
			else if (tag == "camax") tag.value(pm->m_fib.m_camax);
			else if (tag == "beta" ) tag.value(pm->m_fib.m_beta );
			else if (tag == "l0"   ) tag.value(pm->m_fib.m_l0   );
			else if (tag == "refl" ) tag.value(pm->m_fib.m_refl );
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());

		// mark tag as read
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------
// Read rigid materials
//
bool FEBioMaterialSection::ParseRigidMaterial(XMLTag &tag, FERigidMaterial *pm)
{
	FEModel& fem = *GetFEModel();

	if (tag == "center_of_mass") { tag.value(pm->m_rc); pm->m_com = 1; return true; }
	return false;
}

//-----------------------------------------------------------------------------
// Parse FEElasticMixture material 
//
bool FEBioMaterialSection::ParseElasticMixture(XMLTag &tag, FEElasticMixture *pm)
{
	const char* sztype = 0;
	const char* szname = 0;

	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "solid")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an elastic material)
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmat);
		if (pme == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// set the solid material pointer
		pm->AddMaterial(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// TODO: assume that the material becomes stable since it is combined with others
		// in a solid mixture.  (This may not necessarily be true.)
		pme->m_unstable = false;

		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse ParseElasticMultigeneration material 
//
bool FEBioMaterialSection::ParseElasticMultigeneration(XMLTag &tag, FEElasticMultigeneration *pm)
{
	const char* sztype = 0;
	const char* szname = 0;
	const char* szid = 0;
	int id;
	
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "solid")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an elastic material)
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
		{
			clog.printbox("INPUT ERROR", "Invalid elastic solid %s in solid mixture material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// get the growth generation and store it in the material ID
		szid = tag.AttributeValue("gen");
		if (szid) {
			id = atoi(szid) - 1;
			if (id < 0)
				throw XMLReader::InvalidAttributeValue(tag, "gen", szid);
		} else {
			throw XMLReader::MissingAttribute(tag, "gen");
		}
		pme->SetID(id);
		
		// set the solid material pointer
		pm->m_pMat.push_back(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// TODO: assume that the material becomes stable since it is combined with others
		// in a solid mixture.  (This may not necessarily be true.)
		pme->m_unstable = false;
		
		// parse the solid
		ParseMaterial(tag, pmat);
		
		return true;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse FERemodelingElasticMaterial material 
//
bool FEBioMaterialSection::ParseRemodelingSolid(XMLTag &tag, FERemodelingElasticMaterial *pm)
{
	const char* sztype = 0;
	const char* szname = 0;
	
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "solid")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an elastic material)
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
		{
			clog.printbox("INPUT ERROR", "Invalid elastic solid %s in remodeling solid material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid material pointer
		pm->m_pBase = pme;
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	else if (tag == "supply")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. a permeability material)
		FESolidSupply* pme = dynamic_cast<FESolidSupply*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid solid supply %s in remodeling solid material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid supply pointer
		pm->m_pSupp = pme;
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	
	return false;
}


//-----------------------------------------------------------------------------
// Parse FEUncoupledElasticMixture material 
//
bool FEBioMaterialSection::ParseUncoupledElasticMixture(XMLTag &tag, FEUncoupledElasticMixture *pm)
{
	const char* sztype = 0;
	const char* szname = 0;

	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "solid")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an uncoupled elastic material)
		FEUncoupledMaterial* pme = dynamic_cast<FEUncoupledMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
		{
			clog.printbox("INPUT ERROR", "Invalid uncoupled elastic solid %s in uncoupled solid mixture material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid material pointer
		pm->AddMaterial(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// TODO: assume that the material becomes stable since it is combined with others
		// in a solid mixture.  (This may not necessarily be true.)
		pme->m_unstable = false;
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse FEBiphasic material 
//
bool FEBioMaterialSection::ParseBiphasicMaterial(XMLTag &tag, FEBiphasic *pm)
{
	// see if we can find a material property with this name
	int nc = pm->FindComponent(tag.Name());
	if (nc == -1) throw XMLReader::InvalidTag(tag);

	// get the material type
	const char* sztype = tag.AttributeValue("type");

	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// create a new material of this type
	FEBioKernel& febio = FEBioKernel::GetInstance();
	FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// assign the new material to the corresponding material property
	if (pm->SetComponent(nc, pmat) == false)
	{
		clog.printbox("INPUT ERROR: Invalid %s definition in biphasic material %s\n", tag.Name(), pm->GetName());
		throw XMLReader::Error();
	}

	// set the new material's name (if defined)
	if (szname) pmat->SetName(szname);
		
	// parse the new material
	ParseMaterial(tag, pmat);
	
	// done
	return true;
}

//-----------------------------------------------------------------------------
// Parse FEBiphasicSolute material 
//
bool FEBioMaterialSection::ParseBiphasicSoluteMaterial(XMLTag &tag, FEBiphasicSolute *pm)
{
	// see if we can find a material property with this name
	int nc = pm->FindComponent(tag.Name());
	if (nc == -1)
	{
		if (tag == "diffusivity")
		{
			// get the material type
			const char* sztype = tag.AttributeValue("type");
		
			// get the material name
			const char* szname = tag.AttributeValue("name", true);
		
			// create a new material of this type
			FEBioKernel& febio = FEBioKernel::GetInstance();
			FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
			if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
			// make sure the base material is a valid material (i.e. a diffusivity material)
			FESoluteDiffusivity* pme = dynamic_cast<FESoluteDiffusivity*>(pmat);
		
			if (pme == 0)
			{
				clog.printbox("INPUT ERROR", "Invalid diffusivity %s in solute material %s\n", szname, pm->GetName());
				throw XMLReader::Error();
			}
		
			// create a new material of type solute
			if (!pm->m_pSolute) {
				FEMaterial* pmm = febio.Create<FEMaterial>("solute", GetFEModel());
				FESolute* pms = dynamic_cast<FESolute*>(pmm);
				pms->SetSoluteID(0);
				pm->m_pSolute = pms;
				// create a solute with ID 0
				FESoluteData* psd = new FESoluteData;
				psd->m_nID = 0;
				strcpy(psd->m_szname, "neutral");
				FEModel::SetSD(psd);
			}
		
			// set the diffusivity pointer
			pm->m_pSolute->m_pDiff = pme;
		
			// set the material's name
			if (szname) pme->SetName(szname);
		
			// set solute ID in diffusivity
			pme->SetSoluteID(0);
		
			// parse the material
			ParseMaterial(tag, pme);
		
			return true;
		}
		else if (tag == "solubility")
		{
			// get the material type
			const char* sztype = tag.AttributeValue("type");
		
			// get the material name
			const char* szname = tag.AttributeValue("name", true);
		
			// create a new material of this type
			FEBioKernel& febio = FEBioKernel::GetInstance();
			FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
			if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
			// make sure the base material is a valid material (i.e. a solubility material)
			FESoluteSolubility* pme = dynamic_cast<FESoluteSolubility*>(pmat);
		
			if (pme == 0)
			{
				clog.printbox("INPUT ERROR", "Invalid solubility %s in solute material %s\n", szname, pm->GetName());
				throw XMLReader::Error();
			}
		
			// create a new material of type solute
			if (!pm->m_pSolute) {
				FEMaterial* pmm = febio.Create<FEMaterial>("solute", GetFEModel());
				FESolute* pms = dynamic_cast<FESolute*>(pmm);
				pms->SetSoluteID(0);
				pm->m_pSolute = pms;
				// create a solute with ID 0
				FESoluteData* psd = new FESoluteData;
				psd->m_nID = 0;
				strcpy(psd->m_szname, "neutral");
				FEModel::SetSD(psd);
			}
		
			// set the solubility pointer
			pm->m_pSolute->m_pSolub = pme;
		
			// set the material's name
			if (szname) pme->SetName(szname);
		
			// set solute ID in solubility
			pme->SetSoluteID(0);
		
			// parse the material
			ParseMaterial(tag, pme);
		
			return true;
		}
		else throw XMLReader::InvalidTag(tag);
	}

	// get the material type
	// TODO: For solute materials, there is no type
	//       but these materials are registered under "solute" 
	//       so we set the type to that
	const char* sztype = 0;
	if (tag == "solute") sztype = "solute";
	else sztype = tag.AttributeValue("type");

	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// create a new material of this type
	FEBioKernel& febio = FEBioKernel::GetInstance();
	FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// assign the new material to the corresponding material property
	if (pm->SetComponent(nc, pmat) == false)
	{
		clog.printbox("INPUT ERROR: Invalid %s definition in material %s\n", tag.Name(), pm->GetName());
		throw XMLReader::Error();
	}

	// TODO: Solute materials need to have their ID set
	if (dynamic_cast<FESolute*>(pmat))
	{
		FESolute* ps = dynamic_cast<FESolute*>(pmat);
		const char* szid = tag.AttributeValue("sol");
		int nid = atoi(szid) - 1;
		if ((nid < 0) || (nid >= MAX_CDOFS)) throw XMLReader::InvalidAttributeValue(tag, "sol", szid);
		ps->SetSoluteID(nid);
	}

	// set the new material's name (if defined)
	if (szname) pmat->SetName(szname);
		
	// parse the new material
	ParseMaterial(tag, pmat);
	
	return true;
}

//-----------------------------------------------------------------------------
// Parse FESolute material 
//
bool FEBioMaterialSection::ParseSoluteMaterial(XMLTag &tag, FESolute *pm)
{
	// see if we can find a material property with this name
	int nc = pm->FindComponent(tag.Name());
	if (nc == -1) throw XMLReader::InvalidTag(tag);

	// get the material type
	const char* sztype = tag.AttributeValue("type");
		
	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// create a new material of this type
	FEBioKernel& febio = FEBioKernel::GetInstance();
	FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
	if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

	// assign the new material to the corresponding material property
	if (pm->SetComponent(nc, pmat) == false)
	{
		clog.printbox("INPUT ERROR: Invalid %s definition in material %s\n", tag.Name(), pm->GetName());
		throw XMLReader::Error();
	}

	// set the new material's name (if defined)
	if (szname) pmat->SetName(szname);
		
	// parse the new material
	ParseMaterial(tag, pmat);
	
	return true;
}


//-----------------------------------------------------------------------------
// Parse FEChemicalReaction material 
//
bool FEBioMaterialSection::ParseReactionMaterial(XMLTag &tag, FEChemicalReaction *pm)
{
	const char* sztype = 0;
	const char* szname = 0;
	const char* szid = 0;
	
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the stoichiometric coefficients and reaction rates
	if (tag == "vR")
	{
		int id, vR;
		szid = tag.AttributeValue("sbm", true);
		if (szid) {
			id = atoi(szid) - 1;
			if (id < 0)
				throw XMLReader::InvalidAttributeValue(tag, "sbm", szid);
			tag.value(vR);
			pm->SetSolidReactantsCoefficients(id, vR);
		}
		szid = tag.AttributeValue("sol", true);
		if (szid) {
			id = atoi(szid) - 1;
			if ((id < 0) || (id >= MAX_CDOFS))
				throw XMLReader::InvalidAttributeValue(tag, "sol", szid);
			tag.value(vR);
			pm->SetSoluteReactantsCoefficients(id, vR);
		}
		
		return true;
	}
	else if (tag == "vP")
	{
		int id, vP;
		szid = tag.AttributeValue("sbm", true);
		if (szid) {
			id = atoi(szid) - 1;
			if (id < 0)
				throw XMLReader::InvalidAttributeValue(tag, "sbm", szid);
			tag.value(vP);
			pm->SetSolidProductsCoefficients(id, vP);
		}
		szid = tag.AttributeValue("sol", true);
		if (szid) {
			id = atoi(szid) - 1;
			if ((id < 0) || (id >= MAX_CDOFS))
				throw XMLReader::InvalidAttributeValue(tag, "sol", szid);
			tag.value(vP);
			pm->SetSoluteProductsCoefficients(id, vP);
		}
		
		return true;
	}
	else if (tag == "forward_rate")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material
		FEReactionRate* pme = dynamic_cast<FEReactionRate*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid reaction rate %s in multiphasic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the reaction rate pointer
		pm->SetForwardReactionRate(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// parse the material
		ParseMaterial(tag, pme);
		
		return true;
	}
	else if (tag == "reverse_rate")
	{
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material
		FEReactionRate* pme = dynamic_cast<FEReactionRate*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid reaction rate %s in multiphasic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the reaction rate pointer
		pm->SetReverseReactionRate(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
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
	else throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse FETriphasic material 
//
bool FEBioMaterialSection::ParseTriphasicMaterial(XMLTag &tag, FETriphasic *pm)
{
	// we handle solutes differently
	if (tag == "solute")
	{
		const char* sztype = "solute";

		// get the (optional) material name
		const char* szname = tag.AttributeValue("name", true);

		// get the solute ID
		const char* szid = tag.AttributeValue("sol");
		int nid = atoi(szid) - 1;

		// create a new material of this type
		FESolute* psol = new FESolute;
		psol->SetSoluteID(nid);

		// add the solute
		pm->AddSolute(psol);

		// set the new material's name (if defined)
		if (szname) psol->SetName(szname);
		
		// parse the new material
		ParseMaterial(tag, psol);
	}
	else
	{
		// get the material type
		const char* sztype = tag.AttributeValue("type");
		
		// get the (optional) material name
		const char* szname = tag.AttributeValue("name", true);

		// see if we can find a material property with this name
		int nc = pm->FindComponent(tag.Name());
		if (nc == -1) throw XMLReader::InvalidTag(tag);

		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// assign the new material to the corresponding material property
		if (pm->SetComponent(nc, pmat) == false)
		{
			clog.printbox("INPUT ERROR: Invalid %s definition in material %s\n", tag.Name(), pm->GetName());
			throw XMLReader::Error();
		}

		// set the new material's name (if defined)
		if (szname) pmat->SetName(szname);
		
		// parse the new material
		ParseMaterial(tag, pmat);
	}
	
	return true;
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
		// get the material type
		const char* sztype = tag.AttributeValue("type");
		
		// get the material name
		const char* szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. a chemical reaction material)
		FEChemicalReaction* pme = dynamic_cast<FEChemicalReaction*>(pmat);
		
		if (pme == 0)
		{
			clog.printbox("INPUT ERROR", "Invalid chemical reaction %s in multiphasic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the chemical reaction pointer
		pm->AddChemicalReaction(pme);
		
		// set the material's name
		if (szname) pme->SetName(szname);
		
		// parse the material
		ParseMaterial(tag, pme);
	}
	else
	{
		// get the material type
		const char* sztype = tag.AttributeValue("type");
		
		// get the (optional) material name
		const char* szname = tag.AttributeValue("name", true);

		// see if we can find a material property with this name
		int nc = pm->FindComponent(tag.Name());
		if (nc == -1) throw XMLReader::InvalidTag(tag);

		// create a new material of this type
		FEBioKernel& febio = FEBioKernel::GetInstance();
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// assign the new material to the corresponding material property
		if (pm->SetComponent(nc, pmat) == false)
		{
			clog.printbox("INPUT ERROR: Invalid %s definition in material %s\n", tag.Name(), pm->GetName());
			throw XMLReader::Error();
		}

		// set the new material's name (if defined)
		if (szname) pmat->SetName(szname);
		
		// parse the new material
		ParseMaterial(tag, pmat);
	}
	
	return true;
}

//-----------------------------------------------------------------------------
// Parse a viscoelastic material
bool FEBioMaterialSection::ParseViscoElasticMaterial(XMLTag &tag, FEViscoElasticMaterial *pm)
{
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "elastic")
	{
		const char* sztype = 0;
		const char* szname = 0;
		
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an elastic material)
		FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
		{
			clog.printbox("INPUT ERROR", "Invalid elastic solid %s in viscoelastic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}

		// set the solid material pointer
		pm->SetBaseMaterial(pme);
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse an uncoupled viscoelastic material
bool FEBioMaterialSection::ParseUncoupledViscoElasticMaterial(XMLTag &tag, FEUncoupledViscoElasticMaterial *pm)
{
	FEBioKernel& febio = FEBioKernel::GetInstance();
	
	// read the solid material
	if (tag == "elastic")
	{
		const char* sztype = 0;
		const char* szname = 0;
		
		// get the material type
		sztype = tag.AttributeValue("type");
		
		// get the material name
		szname = tag.AttributeValue("name", true);
		
		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, GetFEModel());
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);
		
		// make sure the base material is a valid material (i.e. an uncoupled material)
		FEUncoupledMaterial* pme = dynamic_cast<FEUncoupledMaterial*>(pmat);
		
		// don't allow rigid bodies
		if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
		{
			clog.printbox("INPUT ERROR", "Invalid elastic solid %s in uncoupled viscoelastic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid material pointer
		pm->SetBaseMaterial(pme);
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	
	return false;
}
