//
//  FEBiphasicTangentDiagnostic.h
//  FEBio2
//
//  Created by Gerard Ateshian on 8/20/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBio2__FEBiphasicTangentDiagnostic__
#define __FEBio2__FEBiphasicTangentDiagnostic__

#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
//! The FEBiphasicTangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEBiphasicTangentDiagnostic : public FEDiagnostic
{
public:
    enum TD_Scenario {
        TDS_BIPHASIC_UNIAXIAL
    };
    
public:
    FEBiphasicTangentDiagnostic(FEModel& fem);
    virtual ~FEBiphasicTangentDiagnostic(){}
    
    bool Init();
    
    bool Run();
    
protected:
    void BuildUniaxial();
    
    void deriv_residual(matrix& ke);
    
    void print_matrix(matrix& m);
    
public:
    TD_Scenario	m_scn;
    double		m_strain;
    double      m_pressure;
    double      m_dt;
};

#endif /* defined(__FEBio2__FEBiphasicTangentDiagnostic__) */
