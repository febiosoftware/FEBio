#pragma once
#include "FECore/FEMaterialPoint.h"

//-----------------------------------------------------------------------------
// Define a material point that stores the deformation gradient at previous time point.
class FEViscousMaterialPoint : public FEMaterialPoint
{
public:
    FEViscousMaterialPoint(FEMaterialPoint *pt) : FEMaterialPoint(pt) {}
    
    FEMaterialPoint* Copy();
    
    void Init();

    void Update(const FETimeInfo& timeInfo);

    void Serialize(DumpStream& ar);
    
public:
    mat3d VelocityGradient();
    
    mat3ds RateOfDeformation();
    
public:
    mat3d	m_Fp;		//!< deformation gradient at previous time point
	double	m_dt;		//!< time increment \todo this is a temporary construct. Fix this.
};
