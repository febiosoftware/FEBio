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
    void SetSurface(FESurface* ps);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) {}
    
    //! mark the velocity
    void MarkVelocity();
    
    //! set the velocity
    void SetVelocity();
    
    //! initialization
    bool Init();
    
private:
    double			m_w;        //!< angular speed
    vec3d           m_n;        //!< unit vector along axis of rotation
    vec3d           m_p;        //!< point on axis of rotation
    vector<vec3d>   m_r;        //!< nodal radial positions
    
public:
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FERotationalVelocity_hpp */
