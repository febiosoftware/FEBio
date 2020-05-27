//
//  FEBiphasicFSITraction.hpp
//  FEBioFluid
//
//  Created by Jay Shim on 2/3/20.
//  Copyright © 2020 febio.org. All rights reserved.
//

#pragma once
#include <FECore/FESurfaceLoad.h>
#include "FEFluid.h"

//-----------------------------------------------------------------------------
//! This surface load represents the traction applied on the solid at the
//! interface between a fluid and solid in an FSI analysis.
class FEBIOFLUID_API FEBiphasicFSITraction : public FESurfaceLoad
{
public:
    //! constructor
    FEBiphasicFSITraction(FEModel* pfem);
    
    //! calculate pressure stiffness
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
    //! serialize data
    void Serialize(DumpStream& ar) override;
    
    //! initialization
    bool Init() override;
    
private:
    double GetFluidDilatation(FESurfaceMaterialPoint& mp, double alpha);
    mat3ds GetFluidStress(FESurfaceMaterialPoint& mp);
    
protected:
    vector<double>      m_K;        //!< fluid bulk modulus
    vector<double>      m_s;        //!< scale factor
    vector<bool>        m_bself;    //!< flag if fluid pressure is applied on its own FSI mesh
    vector<FEElement*>  m_elem;     //!< list of fluid-FSI elements
    
    // degrees of freedom
    FEDofList    m_dofU, m_dofSU, m_dofW;
    int        m_dofEF;
};
