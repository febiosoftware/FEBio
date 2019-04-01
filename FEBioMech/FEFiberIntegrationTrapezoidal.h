#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Trapezoidal integration scheme for 2D continuous fiber distributions
//
class FEFiberIntegrationTrapezoidal : public FEFiberIntegrationScheme
{
	class Iterator;

public:
    FEFiberIntegrationTrapezoidal(FEModel* pfem);
    ~FEFiberIntegrationTrapezoidal();

	// get iterator	
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp) override;
    
private:
    int             m_nth;  // number of trapezoidal integration points along theta

	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
