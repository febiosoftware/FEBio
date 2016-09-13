//
//  FEFiberIntegrationGauss.h
//
//  Created by Gerard Ateshian on 11/19/13.
//

#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGauss : public FEFiberIntegrationScheme
{
	class Iterator;

	struct GRULE
	{
		int	m_nph;	// number of gauss integration points along phi
		int	m_nth;  // number of trapezoidal integration points along theta
		const double* m_gp;   // gauss points
		const double* m_gw;   // gauss weights
	};

public:
    FEFiberIntegrationGauss(FEModel* pfem);
    ~FEFiberIntegrationGauss();
	
	//! Initialization
	bool Init();
    
	//! Serialization
	void Serialize(DumpStream& ar);

	// get iterator
	virtual FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp);

protected:
	bool InitRule();
    
protected:	// parameters
	GRULE	m_rule;

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
