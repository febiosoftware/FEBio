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
	static double	m_cth[];
	static double	m_sth[];
	static double	m_cph[];
	static double	m_sph[];
	static double	m_w[];
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
