//
//  FEFluidNormalVelocity.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEFluidNormalVelocity_hpp
#define FEFluidNormalVelocity_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidNormalVelocity is a fluid surface that has a normal
//! velocity prescribed on it.  This routine simultaneously prescribes a
//! surface load and nodal prescribed velocities
class FEFluidNormalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalVelocity(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps);
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R);
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! mark the velocity
    void MarkVelocity();
    
    //! set the velocity
    void SetVelocity();
    
    //! initialization
    bool Init();
    
private:
    double			m_velocity;	//!< magnitude of traction load
    FESurfaceMap	m_VC;		//!< traction boundary cards
    vector<double>  m_VN;       //!< nodal scale factors
    vector<vec3d>   m_nu;       //!< nodal normals

public:
    bool            m_bpv;      //!< flag for prescribing nodal values
    
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidNormalVelocity_hpp */
