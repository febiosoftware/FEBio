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
public:
    FEFiberIntegrationGauss(FEModel* pfem) : FEFiberIntegrationScheme(pfem) { m_nph = 5; m_nth = 2*m_nph; }
    ~FEFiberIntegrationGauss() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
	int             m_nph;	// number of gauss integration points along phi
    int             m_nth;  // number of trapezoidal integration points along theta
    vector<double>  m_gp;   // gauss points
    vector<double>  m_gw;   // gauss weights

	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
