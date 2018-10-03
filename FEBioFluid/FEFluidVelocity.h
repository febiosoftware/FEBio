//
//  FEFluidVelocity.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 4/4/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEFluidVelocity_hpp
#define FEFluidVelocity_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEFluidVelocity is a fluid surface that has a velocity
//! prescribed on it.  This routine simultaneously prescribes a
//! surface load and nodal prescribed velocities
class FEFluidVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidVelocity(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! set the velocity
    void Update() override;
    
    //! initialization
    bool Init() override;
    
    //! activate
    void Activate() override;
    
private:
    double			m_scale;	//!< average velocity
    FESurfaceMap	m_VC;		//!< velocity boundary cards
    vector<vec3d>   m_VN;       //!< nodal velocities
    
public:
    bool            m_bpv;      //!< flag for prescribing nodal values
    
    int		m_dofWX;
    int		m_dofWY;
    int		m_dofWZ;
    int		m_dofEF;
    
    DECLARE_FECORE_CLASS();
};

#endif /* FEFluidVelocity_hpp */
