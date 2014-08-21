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
public:
    FEFiberIntegrationTriangle(FEModel* pfem) : FEFiberIntegrationScheme(pfem) { m_nres = 0; }
    ~FEFiberIntegrationTriangle() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
	int             m_nres;	// resolution
    int             m_nint; // number of integration points
	double          m_cth[2000];
	double          m_sth[2000];
	double          m_cph[2000];
	double          m_sph[2000];
	double          m_w[2000];
    bool            m_bfirst;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
