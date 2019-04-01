#pragma once
#include "FEUncoupledMaterial.h"

//-----------------------------------------------------------------------------
// This material implements an active contraction model which can be used
// as a component of an uncoupled solid mixture material.
class FEPrescribedActiveContractionTransIsoUC : public FEUncoupledMaterial
{
public:
    //! constructor
    FEPrescribedActiveContractionTransIsoUC(FEModel* pfem);
    
	//! validation
	bool Validate() override;

	//! serialization
	void Serialize(DumpStream& ar) override;
    
    //! stress
    mat3ds DevStress(FEMaterialPoint& pt) override;
    
    //! tangent
    tens4ds DevTangent(FEMaterialPoint& pt) override;
    
public: // material parameters
    double	m_thd;		// theta angle for fiber orientation (local coordinates system)
    double	m_phd;		// phi angle for fiber orientation (local coordinates system)
	double	m_T0;       // prescribed active stress

private:
    vec3d	m_n0;		// unit vector along fiber direction (local coordinate system)
    
    DECLARE_FECORE_CLASS();
};
