//
//  FEPrescribedActiveContractionIsotropic.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEPrescribedActiveContractionIsotropic__
#define __FEBioMech__FEPrescribedActiveContractionIsotropic__

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of a solid mixture material.
class FEPrescribedActiveContractionIsotropic : public FEElasticMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionIsotropic(FEModel* pfem);
    
    //! stress
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
public:
    double	m_T0;       // prescribed active stress
    
    DECLARE_FECORE_CLASS();
};


#endif /* defined(__FEBioMech__FEPrescribedActiveContractionIsotropic__) */
