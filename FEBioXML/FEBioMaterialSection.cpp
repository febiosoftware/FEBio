#include "stdafx.h"
#include "FEBioMaterialSection.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
void FEBioMaterialSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();

	// Make sure no materials are defined
	if (fem.Materials() != 0) throw FEFEBioImport::DuplicateMaterialSection();

	const char* sztype = 0;
	const char* szname = 0;

	m_nmat = 0;

	FEBioKernel& febio = FEBioKernel::GetInstance();

	++tag;
	do
	{
		// get the material type
		sztype = tag.AttributeValue("type");

		// get the material name
		szname = tag.AttributeValue("name", true);

		// create a new material of this type
		FEMaterial* pmat = febio.Create<FEMaterial>(sztype, &fem);
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// IMPORTANT: depending on the format version number we need to process 
		// rigid bodies differently. For versions <= 0x0100 rigid degrees of 
		// freedom are initially constrained and can be defined in the material
		// section. For versions >= 0x0101 rigid degrees of freedom are free and
		// can be constrained in the Constraints section.
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(pmat);
		if (pm)
		{
			if (m_pim->Version() <= 0x0100)
			{
				// older versions have the rigid degrees of freedom constrained
				for (int i=0; i<6; ++i) pm->m_bc[i] = -1;
			}
			else
			{
				// newer versions have the rigid degrees of freedom unconstrained
				for (int i=0; i<6; ++i) pm->m_bc[i] = 0;
			}
		}

		// add the material
		fem.AddMaterial(pmat);
		++m_nmat;

		// set the material's name
		if (szname) pmat->SetName(szname);

		// set the material's ID
		pmat->SetID(m_nmat);

		ParseMaterial(tag, pmat);

		// read next tag
		++tag;
	}
	while (!tag.isend());

	// assign material pointers for nested materials
	for (int i=0; i<fem.Materials(); ++i)
	{
		FENestedMaterial* pm = dynamic_cast<FENestedMaterial*>(fem.GetMaterial(i));
		// NOTE: The nested version of visco-elasticity is obselete. Until a new
		//		 implementation is available I am using the FENestedMaterial version
		//		 with some tweaks. Only if the m_nBaseMat is not -1 the nested formulation
		//		 is used. If m_nBaseMat == -1 we assume the new formulation and assume that
		//		 the m_pBase is already set.
		if (pm)
		{
			if (pm->m_nBaseMat == -1)
			{
				if (pm->m_pBase == 0) clog.printbox("INPUT ERROR", "base material for material %d is not defined\n", i+1);
			}
			else
			{
				// get the ID of the base material
				// note that m_nBaseMat is a one-based variable!
				int nbase = pm->m_nBaseMat - 1;

				// make sure the base ID is valid
				if ((nbase < 0) || (nbase >= fem.Materials()))
				{
					clog.printbox("INPUT ERROR", "Invalid base material ID for material %d\n", i+1);
					throw XMLReader::Error();
				}

				// make sure the base material is a valid material (i.e. an elastic material)
				FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(fem.GetMaterial(nbase));

				// don't allow rigid bodies
				if ((pme == 0) || (dynamic_cast<FERigidMaterial*>(pme)))
				{
					clog.printbox("INPUT ERROR", "Invalid base material for material %d\n", i+1);
					throw XMLReader::Error();
				}

				// set the base material pointer
				pm->m_pBase = pme;
			}
		}
	}
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
			
			// triphasic material parameters
			if (!bfound && dynamic_cast<FETriphasic*>(pmat)) bfound = ParseTriphasicMaterial(tag, dynamic_cast<FETriphasic*>(pmat));

			// multiphasic material parameters
			if (!bfound && dynamic_cast<FEMultiphasic*>(pmat)) bfound = ParseMultiphasicMaterial(tag, dynamic_cast<FEMultiphasic*>(pmat));

			// viscoelastic materials
			if (!bfound && dynamic_cast<FEViscoElasticMaterial*>(pmat)) bfound = ParseViscoElasticMaterial(tag, dynamic_cast<FEViscoElasticMaterial*>(pmat));

			// uncoupled viscoelastic materials
			if (!bfound && dynamic_cast<FEUncoupledViscoElasticMaterial*>(pmat)) bfound = ParseUncoupledViscoElasticMaterial(tag, dynamic_cast<FEUncoupledViscoElasticMaterial*>(pmat));
			
			// nested materials
			if (!bfound && dynamic_cast<FENestedMaterial*>(pmat)) bfound = ParseNestedMaterial(tag, dynamic_cast<FENestedMaterial*>(pmat));
			
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
		const char* szt = tag.AttributeValue("type");
		if (strcmp(szt, "local") == 0)
		{
			FELocalMap* pmap = new FELocalMap(mesh);
			pm->m_pmap = pmap;

			int n[3] = {0};
			tag.value(n, 3);
			if ((n[0] == 0) && (n[1] == 0) && (n[2] == 0)) { n[0] = 1; n[1] = 2; n[2] = 4; }

			pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
		}
		else if (strcmp(szt, "vector") == 0)
		{
			FEVectorMap* pmap = new FEVectorMap();
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
		else if (strcmp(szt, "user") == 0)
		{
			// material axis are read in from the ElementData section
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

		return true;
	}
	else if (tag == "fiber")
	{
		const char* szt = tag.AttributeValue("type");
		if (strcmp(szt, "local") == 0)
		{
			FELocalMap* pmap = new FELocalMap(mesh);
			pm->m_pmap = pmap;
			
			int n[3] = {0};
			tag.value(n, 2);
			
			if ((n[0]==0)&&(n[1]==0)&&(n[2]==0)) { n[0] = 1; n[1] = 2; n[2] = 4; }
			if (n[2] == 0) n[2] = n[1];
			
			pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
		}
		else if (strcmp(szt, "spherical") == 0)
		{
			FESphericalMap* pmap = new FESphericalMap(mesh);
			pm->m_pmap = pmap;
			
			vec3d c;
			tag.value(c);
			
			pmap->SetSphereCenter(c);
		}
		else if (strcmp(szt, "cylindrical") == 0)
		{
			FECylindricalMap* pmap = new FECylindricalMap(mesh);
			pm->m_pmap = pmap;
			
			vec3d a(0,0,1), c(0,0,0), r(1,0,0);
			if (tag.isleaf()) throw XMLReader::InvalidValue(tag);

			++tag;
			do
			{
				if      (tag == "center") tag.value(c);
				else if (tag == "axis"  ) tag.value(a);
				else if (tag == "vector") tag.value(r);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			pmap->SetCylinderCenter(c);
			pmap->SetCylinderAxis(a);
			pmap->SetCylinderRef(r);
		}
		else if (strcmp(szt, "vector") == 0)
		{
			FEVectorMap* pmap = new FEVectorMap;
			pm->m_pmap = pmap;
			
			vec3d a, d;
			tag.value(a);
			a.unit();
			
			d = vec3d(1,0,0);
			if (a*d > .999) d = vec3d(0,1,0);
			
			pmap->SetVectors(a, d);
		}
		else if (strcmp(szt, "user") == 0)
		{
			// fibers are read in in the ElementData section
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		
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
		const char* szt = tag.AttributeValue("type");
		if (strcmp(szt, "local") == 0)
		{
			FELocalMap* pmap = new FELocalMap(mesh);
			pm->m_pmap = pmap;

			int n[3] = {0};
			tag.value(n, 2);

			if ((n[0]==0)&&(n[1]==0)&&(n[2]==0)) { n[0] = 1; n[1] = 2; n[2] = 4; }
			if (n[2] == 0) n[2] = n[1];

			pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
		}
		else if (strcmp(szt, "spherical") == 0)
		{
			FESphericalMap* pmap = new FESphericalMap(mesh);
			pm->m_pmap = pmap;

			vec3d c;
			tag.value(c);

			pmap->SetSphereCenter(c);
		}
		else if (strcmp(szt, "cylindrical") == 0)
		{
			FECylindricalMap* pmap = new FECylindricalMap(mesh);
			pm->m_pmap = pmap;
			
			vec3d a(0,0,1), c(0,0,0), r(1,0,0);
			if (tag.isleaf()) throw XMLReader::InvalidValue(tag);

			++tag;
			do
			{
				if      (tag == "center") tag.value(c);
				else if (tag == "axis"  ) tag.value(a);
				else if (tag == "vector") tag.value(r);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			pmap->SetCylinderCenter(c);
			pmap->SetCylinderAxis(a);
			pmap->SetCylinderRef(r);
		}
		else if (strcmp(szt, "vector") == 0)
		{
			FEVectorMap* pmap = new FEVectorMap;
			pm->m_pmap = pmap;

			vec3d a, d;
			tag.value(a);
			a.unit();

			d = vec3d(1,0,0);
			if (a*d > .999) d = vec3d(0,1,0);

			pmap->SetVectors(a, d);
		}
		else if (strcmp(szt, "user") == 0)
		{
			// fibers are read in in the ElementData section
		}
		else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

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
			if (tag == "ca0") tag.value(pm->m_fib.m_ca0);
			else if (tag == "beta") tag.value(pm->m_fib.m_beta);
			else if (tag == "l0") tag.value(pm->m_fib.m_l0);
			else if (tag == "refl") tag.value(pm->m_fib.m_refl);
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
	else if (tag == "parent_id") { tag.value(pm->m_pmid); return true; }
	else if (m_pim->Version() <= 0x0100)
	{
		FEAnalysisStep* pStep = GetStep();

		// The following tags are only allowed in older version of FEBio
		// Newer versions defined the rigid body constraints in the Constraints section
		if (strncmp(tag.Name(), "trans_", 6) == 0)
		{
			const char* szt = tag.AttributeValue("type");

			int bc = -1;
			if      (tag.Name()[6] == 'x') bc = 0;
			else if (tag.Name()[6] == 'y') bc = 1;
			else if (tag.Name()[6] == 'z') bc = 2;
			assert(bc >= 0);

			if      (strcmp(szt, "free"      ) == 0) pm->m_bc[bc] =  0;
			else if (strcmp(szt, "fixed"     ) == 0) pm->m_bc[bc] = -1;
			else if (strcmp(szt, "prescribed") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc)-1;

				pm->m_bc[bc] = lc+1;
				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement;
				pDC->id = m_nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RDC.size()-1;
					FERigidBodyDisplacement* pDC = fem.m_RDC[n];
					pStep->AddBoundaryCondition(pDC);
					pDC->Deactivate();
				}
			}
			else if (strcmp(szt, "force") == 0)
			{
				const char* szlc = tag.AttributeValue("lc");
				int lc = atoi(szlc)-1;

				pm->m_bc[bc] = 0;
				FERigidBodyForce* pFC = new FERigidBodyForce;
				pFC->id = m_nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
			return true;
		}
		else if (strncmp(tag.Name(), "rot_", 4) == 0)
		{
			const char* szt = tag.AttributeValue("type");

			int bc = -1;
			if      (tag.Name()[4] == 'x') bc = 3;
			else if (tag.Name()[4] == 'y') bc = 4;
			else if (tag.Name()[4] == 'z') bc = 5;
			assert(bc >= 0);

			if      (strcmp(szt, "free"      ) == 0) pm->m_bc[bc] =  0;
			else if (strcmp(szt, "fixed"     ) == 0) pm->m_bc[bc] = -1;
			else if (strcmp(szt, "prescribed") == 0)
			{
				const char* szlc = tag.AttributeValue("lc", true);
				int lc = atoi(szlc)-1;

				pm->m_bc[bc] = lc+1;
				FERigidBodyDisplacement* pDC = new FERigidBodyDisplacement;
				pDC->id = m_nmat;
				pDC->bc = bc;
				pDC->lc = lc;
				tag.value(pDC->sf);
				fem.m_RDC.push_back(pDC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RDC.size()-1;
					FERigidBodyDisplacement* pDC = fem.m_RDC[n];
					pStep->AddBoundaryCondition(pDC);
					pDC->Deactivate();
				}
			}
			else if (strcmp(szt, "force") == 0)
			{
				const char* szlc = tag.AttributeValue("lc", true);
				int lc = atoi(szlc)-1;

				pm->m_bc[bc] = 0;
				FERigidBodyForce* pFC = new FERigidBodyForce;
				pFC->id = m_nmat;
				pFC->bc = bc;
				pFC->lc = lc;
				tag.value(pFC->sf);
				fem.m_RFC.push_back(pFC);

				// add this boundary condition to the current step
				if (m_pim->m_nsteps > 0)
				{
					int n = fem.m_RFC.size()-1;
					FERigidBodyForce* pFC = fem.m_RFC[n];
					pStep->AddBoundaryCondition(pFC);
					pFC->Deactivate();
				}
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
			return true;
		}
	}
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
		pm->m_pMat.push_back(pme);
		
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
	if (nc == -1) throw XMLReader::InvalidTag(tag);

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
// Parse FETriphasic material 
//
bool FEBioMaterialSection::ParseTriphasicMaterial(XMLTag &tag, FETriphasic *pm)
{
	// get the material type
	const char* sztype = tag.AttributeValue("type");
		
	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// get the material ID
	// TODO: If this is a triphasic material, this is only used for the solutes.
	//       Perhaps I can encode whether the ID is necessary in the FEMultiMaterial class.
	//       For now I have to define it as an optional argument.
	int nid = 0;
	const char* szid = tag.AttributeValue("sol", true);
	if (szid) nid = atoi(szid) - 1;

	// see if we can find a material property with this name
	int nc = pm->FindComponent(tag.Name(), nid);
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
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse FETriphasic material 
//
bool FEBioMaterialSection::ParseMultiphasicMaterial(XMLTag &tag, FEMultiphasic *pm)
{
	// get the material type
	const char* sztype = tag.AttributeValue("type");
		
	// get the material name
	const char* szname = tag.AttributeValue("name", true);

	// get the material ID
	// TODO: If this is a triphasic material, this is only used for the solutes.
	//       Perhaps I can encode whether the ID is necessary in the FEMultiMaterial class.
	//       For now I have to define it as an optional argument.
	int nid = 0;
	const char* szid = tag.AttributeValue("sol", true);
	if (szid) nid = atoi(szid) - 1;

	// see if we can find a material property with this name
	int nc = pm->FindComponent(tag.Name(), nid);
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
	
	return false;
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
		pm->m_pBase = pme;
		
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
		pm->m_pBase = pme;
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}
	
	return false;
}

//-----------------------------------------------------------------------------
// Parse a nested material
// NOTE: nested materials are obselete, but we are still using them in the
// interest of backward compatibility. At this point, they are somewhat of a
// hybrid between the old and new format.
bool FEBioMaterialSection::ParseNestedMaterial(XMLTag &tag, FENestedMaterial *pm)
{
	FEBioKernel& febio = FEBioKernel::GetInstance();

	// Make sure the m_nBaseMat is -1
	if (pm->m_nBaseMat != -1) return false;
	
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
			clog.printbox("INPUT ERROR", "Invalid elastic solid %s in biphasic material %s\n", szname, pm->GetName());
			throw XMLReader::Error();
		}
		
		// set the solid material pointer
		pm->m_pBase = pme;
		
		// parse the solid
		ParseMaterial(tag, pme);
		
		return true;
	}

	return false;
}
