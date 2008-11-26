#include "stdafx.h"
#include "FEMaterialFactory.h"

// This is the one and only instance of the material factory
// You can retrieve this pointer by calling FEMaterialFactory::GetInstance()
FEMaterialFactory*	FEMaterialFactory::m_pMF = 0;

//-----------------------------------------------------------------------------
//! This function can be used to retrieve the one and only instance of the
//! FEMaterialFactory class.
FEMaterialFactory* FEMaterialFactory::GetInstance()
{
	if (m_pMF == 0) m_pMF = new FEMaterialFactory;
	assert(m_pMF);
	return m_pMF;
}

//-----------------------------------------------------------------------------
//! This function is used to register a name with a material. Or to be more 
//! precise, to associate a name with the function that will create this
//! material. This function is usually a static member function of the material
//! class. Users should usually not want to call this function directly.
//! Instead they should add the DECLARE_REGISTERED and REGISTER_MATERIAL macros.
void FEMaterialFactory::RegisterMaterial(const char* sz, FEMaterialFactory::MAT_CREATE_FUNC fnc)
{
	// create a REGISTER_MATERIAL object
	REGISTERED_MATERIAL rm;

	// get the instance
	FEMaterialFactory* pMF = GetInstance();

	// get the material list
	list<REGISTERED_MATERIAL>& ml = pMF->m_ml;

	// search the list to make sure that the material name has not been registered yet
	if (ml.size()>0)
	{
		list<REGISTERED_MATERIAL>::iterator it;
		for (it = ml.begin(); it != ml.end(); ++it)
			if (strcmp(it->m_szname, sz) == 0)
			{
				fprintf(stderr, "FATAL ERROR: a material has already been declared with this name.\n\n");
				throw FEMaterialFactory::Exception();
			}
	}

	// make sure the name is not too long
	if (strlen(sz)>=MAX_MATERIAL_NAME)
	{
		fprintf(stderr, "FATAL ERROR: registered material name is too long.\n\n");
		throw FEMaterialFactory::Exception();
	}

	// set the parameters
	rm.m_fnc = fnc;
	strcpy(rm.m_szname, sz);

	// add the registration to the list
	ml.push_back(rm);
}

//-----------------------------------------------------------------------------
//! This function is used by the framework to create a specific material
//! from its registered name.
//! \param szmat the registered name of the material
FEMaterial* FEMaterialFactory::CreateMaterial(const char *szmat)
{
	// get the instance
	FEMaterialFactory* pMF = GetInstance();

	// get the material list
	list<REGISTERED_MATERIAL>& ml = pMF->m_ml;

	// find the name
	if (ml.size() > 0)
	{
		list<REGISTERED_MATERIAL>::iterator it;
		for (it = ml.begin(); it != ml.end(); ++it)
		{
			if (strcmp(it->m_szname, szmat) == 0) return it->m_fnc();
		}
	}

	// if we get here then the material was not found
	return 0;
}
