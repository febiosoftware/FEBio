//
//  FEMultiphasicTangentDiagnostic.h
//  FEBio2
//
//  Created by Gerard Ateshian on 8/21/15.
//  Copyright (c) 2015 febio.org. All rights reserved.
//

#ifndef __FEBio2__FEMultiphasicTangentDiagnostic__
#define __FEBio2__FEMultiphasicTangentDiagnostic__

#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
//! The FEMultiphasicTangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEMultiphasicTangentDiagnostic : public FEDiagnostic
{
public:
    enum TD_Scenario {
        TDS_MULTIPHASIC_UNIAXIAL
    };
    
public:
    FEMultiphasicTangentDiagnostic(FEModel& fem);
    virtual ~FEMultiphasicTangentDiagnostic(){}
    
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
    double      m_concentration;
    double      m_dt;
};

#endif /* defined(__FEBio2__FEMultiphasicTangentDiagnostic__) */
