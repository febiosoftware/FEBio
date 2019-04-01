#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of a solid mixture material.
class FEPrescribedActiveContractionIsotropicUC : public FEUncoupledMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionIsotropicUC(FEModel* pfem);
    
    //! stress
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! tangent
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
public:
    double	m_T0;       // prescribed active stress
    
    DECLARE_FECORE_CLASS();
};
