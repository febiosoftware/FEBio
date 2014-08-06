//
//  FEFiberIntegrationTriangleUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 8/5/14.
//  Copyright (c) 2014 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEFiberIntegrationTriangleUC__
#define __FEBioMech__FEFiberIntegrationTriangleUC__

#include "FEFiberIntegrationSchemeUC.h"

//----------------------------------------------------------------------------------
// Geodesic dome integration scheme for continuous fiber distributions
//
class FEFiberIntegrationTriangleUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationTriangleUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nres = 0; }
    ~FEFiberIntegrationTriangleUC() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds DevStress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds DevTangent(FEMaterialPoint& mp);
    
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

#endif /* defined(__FEBioMech__FEFiberIntegrationTriangleUC__) */
