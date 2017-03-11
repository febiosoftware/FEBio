//
//  FEBackFlowStabilization.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 3/2/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEBackFlowStabilization_hpp
#define FEBackFlowStabilization_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>

//-----------------------------------------------------------------------------
//! Backflow stabilization prescribes a normal traction that opposes
//! backflow on a boundary surface.
class FEBackFlowStabilization : public FESurfaceLoad
{
public:
    //! constructor
    FEBackFlowStabilization(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) {}
    
    //! serialize data
    void Serialize(DumpStream& ar);
    
    //! mark the dilatation
    void MarkDilatation();
    
    //! set the dilatation
    void SetDilatation();
    
    //! initialization
    bool Init();
    
protected:
    //! calculate stiffness for an element
    void ElementStiffness(FESurfaceElement& el, matrix& ke);
    
    //! Calculates the force for an element
    void ElementForce(FESurfaceElement& el, vector<double>& fe);
    
protected:
    double			m_beta;     //!< backflow stabilization coefficient
    double          m_rho;      //!< fluid density
    bool            m_bout;     //!< outflow flag
    double          m_dir;      //!< flow direction
    double          m_k;        //!< fluid bulk modulus
    vector<vec3d>   m_nu;       //!< nodal normals
    
    // degrees of freedom
    int		m_dofVX;
    int		m_dofVY;
    int		m_dofVZ;
    int		m_dofE;
    
    DECLARE_PARAMETER_LIST();
};


#endif /* FEBackFlowStabilization_hpp */
