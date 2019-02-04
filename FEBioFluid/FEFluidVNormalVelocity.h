//
//  FEFluidVNormalVelocity.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/2/19.
//  Copyright © 2019 febio.org. All rights reserved.
//

#ifndef FEFluidVNormalVelocity_hpp
#define FEFluidVNormalVelocity_hpp

#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FEFluidVNormalVelocity is a fluid surface that has a normal
//! velocity prescribed on it.  This routine prescribes
//! nodal velocities with option of a parabolic profile
class FEBIOFLUID_API FEFluidVNormalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidVNormalVelocity(FEModel* pfem);
    
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
    
    //! parabolic velocity profile
    bool SetParabolicVelocity();
    
private:
    double          m_velocity;    //!< average velocity
    FESurfaceMap    m_VC;        //!< velocity boundary cards
    vector<double>  m_VN;       //!< nodal scale factors
    vector<vec3d>   m_nu;       //!< nodal normals
    bool            m_bpar;     //!< flag for parabolic velocity
    
public:
    int        m_dofWX;
    int        m_dofWY;
    int        m_dofWZ;
    
    DECLARE_FECORE_CLASS();
};

#endif /* FEFluidVNormalVelocity_hpp */
