#pragma once
#include "FEMultiphasic.h"

//-----------------------------------------------------------------------------
//! Standard multiphasic material.

class FEMultiphasicStandard : public FEMultiphasic
{
public:
	//! constructor
	FEMultiphasicStandard(FEModel* pfem);
    
    //! returns a pointer to a new material point object
	FEMaterialPoint* CreateMaterialPointData();
	
    //! Update solid bound molecules
    void UpdateSolidBoundMolecules(FEMaterialPoint& mp);
};
