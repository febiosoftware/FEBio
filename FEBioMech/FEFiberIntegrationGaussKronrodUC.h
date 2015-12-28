//
//  FEFiberIntegrationGaussKronrodUC.h
//
//

#pragma once
#include "FEFiberIntegrationSchemeUC.h"

//----------------------------------------------------------------------------------
// Gauss integration scheme for continuous fiber distributions with uncoupled strain energy
//
class FEFiberIntegrationGaussKronrodUC : public FEFiberIntegrationSchemeUC
{
public:
    FEFiberIntegrationGaussKronrodUC(FEModel* pfem) : FEFiberIntegrationSchemeUC(pfem) { m_nph = 7; m_nth = 31; }
    ~FEFiberIntegrationGaussKronrodUC() {}
    
    //! Initialization
    bool Init();
    
    //! Cauchy stress
    mat3ds DevStress(FEMaterialPoint& mp);
    
    // Spatial tangent
    tens4ds DevTangent(FEMaterialPoint& mp);
    
    //! Strain energy density
    double DevStrainEnergyDensity(FEMaterialPoint& mp);
    
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
