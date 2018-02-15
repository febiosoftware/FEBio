//
//  FETangentialDamping.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/1/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FETangentialDamping_hpp
#define FETangentialDamping_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! Tangential damping prescribes a shear traction that opposes tangential
//! fluid velocity on a boundary surface.  This can help stabilize inflow
//! conditions.
class FETangentialDamping : public FESurfaceLoad
{
public:
    //! constructor
    FETangentialDamping(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override;
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const double alpha);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const double alpha);
    
protected:
    double			m_eps;      //!< damping coefficient (penalty)
    
    // degrees of freedom
    int		m_dofWX, m_dofWY, m_dofWZ;
    int		m_dofWXP, m_dofWYP, m_dofWZP;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FETangentialDamping_hpp */
