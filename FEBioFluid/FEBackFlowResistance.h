//
//  FEBackFlowResistance.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 8/8/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEBackFlowResistance_hpp
#define FEBackFlowResistance_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! FEBackflowResistance is a fluid surface that has a normal
//! pressure that resists backflow.
//!
class FEBackFlowResistance : public FESurfaceLoad
{
public:
    //! constructor
    FEBackFlowResistance(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override { m_alpha = tp.alpha; m_alphaf = tp.alphaf; }
    
    //! set the dilatation
    void Update() override;
    
    //! evaluate normal velocities
    void NormalVelocities();
    
    //! evaluate flow rate
    double FlowRate();
    
    //! initialize
    bool Init() override;
    
    //! activate
    void Activate() override;
    
private:
    double          m_beta;     //!< flow resistance
    double          m_R;        //!< flow resistance
    double          m_k;        //!< fluid bulk modulus
    double          m_rho;      //!< fluid density
    double          m_alpha;
    double          m_alphaf;
    double          m_p0;       //!< fluid pressure offset
    vector<double>  m_vn;       //!< nodal normal velocity
    
    int        m_dofWX, m_dofWY, m_dofWZ;
    int        m_dofWXP, m_dofWYP, m_dofWZP;
    int        m_dofEF;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEBackFlowResistance_hpp */
