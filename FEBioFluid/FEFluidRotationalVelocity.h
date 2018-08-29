//
//  FERotationalVelocity.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 6/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FERotationalVelocity_hpp
#define FERotationalVelocity_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidRotationalVelocity is a fluid surface that has a rotational
//! velocity prescribed on it.  This routine prescribes nodal velocities

class FEFluidRotationalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidRotationalVelocity(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override {}
    
    //! set the velocity
    void Update() override;
    
    //! initialization
    bool Init() override;
    
    //! activate
    void Activate() override;
    
private:
    double			m_w;        //!< angular speed
    vec3d           m_n;        //!< unit vector along axis of rotation
    vec3d           m_p;        //!< point on axis of rotation
    vector<vec3d>   m_r;        //!< nodal radial positions
    
public:
    int		m_dofWX;
    int		m_dofWY;
    int		m_dofWZ;
    int		m_dofEF;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FERotationalVelocity_hpp */
