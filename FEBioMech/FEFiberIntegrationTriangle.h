//
//  FEFiberIntegrationTriangle.h
//  FEBioXCode4
//
//

#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Geodesic dome integration scheme for continuous fiber distributions
//
class FEFiberIntegrationTriangle : public FEFiberIntegrationScheme
{
	class Iterator;

public:
    FEFiberIntegrationTriangle(FEModel* pfem);
    ~FEFiberIntegrationTriangle();
	
	//! Initialization
	bool Init();
    
	// serialization
	void Serialize(DumpStream& ar);

	// create iterator
	FEFiberIntegrationSchemeIterator* GetIterator(FEMaterialPoint* mp);

protected:
	void InitIntegrationRule();
    
public: // parameters
	int             m_nres;	// resolution

protected:
    int             m_nint; // number of integration points
	double          m_cth[2000];
	double          m_sth[2000];
	double          m_cph[2000];
	double          m_sph[2000];
	double          m_w[2000];
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
