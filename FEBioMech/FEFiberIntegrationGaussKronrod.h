#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGaussKronrod : public FEFiberIntegrationScheme
{
public:
	class Iterator;

	struct GKRULE
	{
		int	m_nph; // number of gauss integration points along phi
		int	m_nth; // number of trapezoidal integration points along theta
		const double*	m_gp; // gauss points
		const double*	m_gw; // gauss weights
	};

public:
    FEFiberIntegrationGaussKronrod(FEModel* pfem);
    ~FEFiberIntegrationGaussKronrod();
	
	//! Initialization
	bool Init() override;
    
	// Serialization
	void Serialize(DumpStream& ar) override;

	// get the iterator
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp) override;

protected:
	bool InitRule();
    
protected: // parameters
	GKRULE	m_rule;
    
	// declare the parameter list
	DECLARE_FECORE_CLASS();
};
