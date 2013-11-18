#pragma once
#include "febio.h"
#include <list>
using namespace std;

class FEMaterial;

typedef FEBioFactory_T<FEMaterial> FEMaterialFactory;

//---------------------------------------------------------------------------------
//! This class helps with the registration of a material with the FEMaterialFactory
// NOTE that I'm not using the REGISTER_FEBIO_CLASS for registering materials since
// this macro assumes that the constructor of the derived class takes a pointer to
// a FEModel class, which is not the case for FEMaterial's. I haven't decided if I
// want to make that the default for all classes yet.
template <typename M> class FERegisterMaterial_T : public FEMaterialFactory
{
public:
	FERegisterMaterial_T(const char* sz) : FEMaterialFactory(sz)
	{
		FEBioKernel& febio = FEBioKernel::GetInstance();
		febio.RegisterClass(this);
	}

	M* Create(FEModel* pfem) { return new M(pfem); }
	bool IsType(FEMaterial* pm) { return (dynamic_cast<M*>(pm) != 0); }
};

// The DECLARE_REGISTERED macro sets up the mechanism to do the material registration
// (this doesn't do anything anymore in FEBio2)
#define DECLARE_REGISTERED(theClass)

// the REGISTER_MATERIAL does the actual material registration
#define REGISTER_MATERIAL(theClass, theName) \
	static FERegisterMaterial_T<theClass> _##theClass##_rm(theName);
