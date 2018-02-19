//
//  FEFluidBCFormula.hpp
//  FEBioFluid
//
//  Created by Gerard Ateshian on 2/18/18.
//  Copyright © 2018 febio.org. All rights reserved.
//

#ifndef FEFluidBCFormula_hpp
#define FEFluidBCFormula_hpp

#include "FECore/FESurfaceLoad.h"
#include <FECore/FESurfaceMap.h>
#include <FECore/MathParser.h>

//-----------------------------------------------------------------------------
//! FEFluidBCFormula is a fluid surface on which fluid nodal boundary
//! conditions may be prescribed using a formula.
class FEFluidBCFormula : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidBCFormula(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(const FETimeInfo& tp, FESolver* psolver) override {}
    
    //! calculate residual
    void Residual(const FETimeInfo& tp, FEGlobalVector& R) override {}
    
    //! mark the boundary condition
    void MarkBC();
    
    //! set the boundary condition
    void SetBC();
    
    //! initialization
    bool Init() override;
    
private:
    double      m_s;        //! scale factor
    char        m_bc[3];    //! boundary conditions
    char        m_sz[256];  //! formula
    int         m_dof;      //! boundary condition DOF

public:
    int         m_dofWX;
    int         m_dofWY;
    int         m_dofWZ;
    int         m_dofEF;
    
    DECLARE_PARAMETER_LIST();
};

#endif /* FEFluidBCFormula_hpp */
