//
//  FEPrescribedActiveContractionTransIso.h
//  FEBioMech
//
//  Created by Gerard Ateshian on 3/13/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBioMech__FEPrescribedActiveContractionTransIso__
#define __FEBioMech__FEPrescribedActiveContractionTransIso__

#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of a solid mixture material.
class FEPrescribedActiveContractionTransIso : public FEElasticMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionTransIso(FEModel* pfem);
    
    //! Validation
    bool Validate();

	//! serialization
	void Serialize(DumpStream& ar);
    
    //! stress
    mat3ds Stress(FEMaterialPoint& pt);
    
    //! tangent
    tens4ds Tangent(FEMaterialPoint& pt);
    
public:
    double	m_T0;       // prescribed active stress
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
    
    DECLARE_PARAMETER_LIST();
};

#endif /* defined(__FEBioMech__FEPrescribedActiveContractionTransIso__) */
