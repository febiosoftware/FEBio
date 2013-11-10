#pragma once
#include "FEUncoupledMaterial.h"
#include "FEFiberMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for transversely isotropic materials.

//! This class was created to simplify the implementation of the TransIso Mooney-
//! Rivlin and Veronda-Westmann models.
//! This class only stores some data

//! \todo Don't derive it from FEUncoupledMaterial. Materials that are both
//!       incompressible and have a fiber distribution should derive from 
//!       FEUncoupledMaterial and from a fiber material class. Or perhaps they should
//!       have a fiber (or material axis) class as a member.

class FETransverselyIsotropic : public FEUncoupledMaterial
{
public:
	//! constructor
	FETransverselyIsotropic() {}

	//! Initialization
	void Init();

	//! serialize material data
	void Serialize(DumpFile& ar);

public:
	//! return the number of properties
	int Properties();

	//! return a pointer to the property
	FEMaterial* GetProperty(int n);

	//! find a material property index ( returns <0 for error)
	int FindPropertyIndex(const char* szname);

	//! set a material property (returns false on error)
	bool SetProperty(int i, FEMaterial* pm);

public:
	FEFiberMaterial	m_fib;

	// declare parameter list
	DECLARE_PARAMETER_LIST();
};
