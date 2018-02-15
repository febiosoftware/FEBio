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
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! mark the velocity
    void MarkVelocity();
    
    //! set the velocity
    void SetVelocity();
    
    //! initialization
    bool Init() override;
    
    //! parabolic velocity profile
    bool SetParabolicVelocity();
    
private:
    double			m_velocity;	//!< average velocity
    FESurfaceMap	m_VC;		//!< velocity boundary cards
    vector<double>  m_VN;       //!< nodal scale factors
    vector<vec3d>   m_nu;       //!< nodal normals
    bool            m_bpar;     //!< flag for parabolic velocity

public:
    bool            m_bpv;      //!< flag for prescribing nodal values
    
    int		m_dofWX;
    int		m_dofWY;
    int		m_dofWZ;
    int		m_dofEF;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidNormalVelocity_hpp */
