#pragma once
#include "FECore/FEMaterial.h"

//-----------------------------------------------------------------------------
//! Base class for the elastic part of the fluid response.
//! These materials provide a relation between fluid pressure and J = det F.
//!
class FEElasticFluid : public FEMaterial
{
public:
    FEElasticFluid(FEModel* pfem) : FEMaterial(pfem) {}
	virtual ~FEElasticFluid() {}
    
	//! fluid pressure
    virtual double Pressure(FEMaterialPoint& pt) = 0;
    
	//! tangent of fluid pressure with respect to strain J
    virtual double Tangent_Pressure_Strain(FEMaterialPoint& mp) = 0;
    
    //! 2nd derivative of fluid pressure with respect to strain J
    virtual double Tangent_Pressure_Strain_Strain(FEMaterialPoint& mp) = 0;
    
    //! bulk modulus
    double BulkModulus(FEMaterialPoint& mp);
};
