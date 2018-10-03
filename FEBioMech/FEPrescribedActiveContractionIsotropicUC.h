//
//  FEPrescribedActiveContractionIsotropicUC.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEPrescribedActiveContractionIsotropicUC__
#define __FEBioMech__FEPrescribedActiveContractionIsotropicUC__

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

#endif /* defined(__FEBioMech__FEPrescribedActiveContractionIsotropicUC__) */
