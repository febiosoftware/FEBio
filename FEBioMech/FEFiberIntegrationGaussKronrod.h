//
//  FEFiberIntegrationGaussKronrod.h
//
//

#pragma once
#include "FEFiberIntegrationScheme.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions
//
class FEFiberIntegrationGaussKronrod : public FEFiberIntegrationScheme
{
public:
    FEFiberIntegrationGaussKronrod(FEModel* pfem) : FEFiberIntegrationScheme(pfem) { m_nph = 7; m_nth = 31; }
    ~FEFiberIntegrationGaussKronrod() {}
	
	//! Initialization
	void Init();
    
	//! Cauchy stress
	mat3ds Stress(FEMaterialPoint& mp);
    
	// Spatial tangent
	tens4ds Tangent(FEMaterialPoint& mp);
    
	//! Strain energy density
	double StrainEnergyDensity(FEMaterialPoint& mp);
    
    // Fiber density
    void IntegratedFiberDensity(double& IFD);
    
public:
	int             m_nph;	// number of gauss integration points along phi
    int             m_nth;  // number of trapezoidal integration points along theta
    vector<double>  m_gp;   // gauss points
    vector<double>  m_gw;   // gauss weights
    bool            m_bfirst;
    
	// declare the parameter list
	DECLARE_PARAMETER_LIST();
};
