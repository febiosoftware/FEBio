#pragma once
#include "FEActiveMomentumSupply.h"

//-----------------------------------------------------------------------------
// This class implements a constant active momentum supply

class FEBIOMIX_API FEActiveConstantSupply :	public FEActiveMomentumSupply
{
public:
    //! constructor
    FEActiveConstantSupply(FEModel* pfem);
    
    //! active momentum supply
    vec3d ActiveSupply(FEMaterialPoint& pt) override;
    
    //! Tangent of active momentum supply
    vec3d Tangent_ActiveSupply_Strain(FEMaterialPoint& mp) override;
        
public:
    double	m_asupp;			//!< active momentum supply
    
    // declare parameter list
    DECLARE_FECORE_CLASS();
};
