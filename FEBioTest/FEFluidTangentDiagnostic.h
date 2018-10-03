//
//  FEFluidTangentDiagnostic.hpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 11/17/15.
//  Copyright © 2015 febio.org. All rights reserved.
//

#ifndef FEFluidTangentDiagnostic_h
#define FEFluidTangentDiagnostic_h

#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FEFluidScenario : public FEDiagnosticScenario
{
public:
    FEFluidScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }
    
public:
    double	m_dt;
};

//-----------------------------------------------------------------------------
class FEFluidTangentUniaxial : public FEFluidScenario
{
public:
    FEFluidTangentUniaxial(FEDiagnostic* pdia);
    
    bool Init() override;
    
private:
    double		m_velocity;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEFluidTangentUniaxialSS : public FEFluidScenario
{
public:
    FEFluidTangentUniaxialSS(FEDiagnostic* pdia);
    
    bool Init() override;
    
private:
    double		m_velocity;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! The FEBiphasicTangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEFluidTangentDiagnostic : public FEDiagnostic
{
public:
    FEFluidTangentDiagnostic(FEModel& fem);
    virtual ~FEFluidTangentDiagnostic(){}
    
    bool Init();
    
    bool Run();
    
    FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void deriv_residual(matrix& ke);
    
    void print_matrix(matrix& m);
    
public:
    FEFluidScenario*	m_pscn;
};

#endif /* FEFluidTangentDiagnostic_h */
