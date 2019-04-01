#pragma once
#include "FEElasticMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of a solid mixture material.
class FEPrescribedActiveContractionUniaxial : public FEElasticMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionUniaxial(FEModel* pfem);
    
    //! validation
    bool Validate() override;

	//! serialization
	void Serialize(DumpStream& ar) override;
    
    //! stress
    mat3ds Stress(FEMaterialPoint& pt) override;
    
    //! tangent
    tens4ds Tangent(FEMaterialPoint& pt) override;
    
public:
    double	m_T0;       // prescribed active stress
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)

private:
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
    
    DECLARE_FECORE_CLASS();
};
