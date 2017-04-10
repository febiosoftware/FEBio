//
//  FETangentialFlowStabilization.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FETangentialFlowStabilization_hpp
#define FETangentialFlowStabilization_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! Tangential flow stabilization prescribes a shear traction that opposes
//! tangential fluid velocity on a boundary surface, in the presence of normal
//! flow.  This can help stabilize inflow/outflow conditions.
class FETangentialFlowStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FETangentialFlowStabilization(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver);
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R);
    
    //! serialize data
    void Serialize(DumpStream& ar);
    
    //! Unpack surface element data
    void UnpackLM(FEElement& el, vector<int>& lm);
    
    //! initialization
    bool Init();
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke, const double alpha);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe, const double alpha);
    
protected:
    double			m_beta;     //!< damping coefficient
    double          m_rho;      //!< fluid density
    
    // degrees of freedom
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    
    DECLARE_PARAMETER_LIST();
};


#endif /* FETangentialFlowStabilization_hpp */
