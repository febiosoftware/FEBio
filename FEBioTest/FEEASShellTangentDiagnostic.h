//
//  FEShellTangentDiagnostic.hpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 11/27/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEShellTangentDiagnostic_hpp
#define FEShellTangentDiagnostic_hpp

#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FEEASShellTangentUnloaded : public FEDiagnosticScenario
{
public:
    FEEASShellTangentUnloaded(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_strain = 0.0; }
    
    bool Init() override;
    
private:
    double    m_strain;
    
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
//! The FETangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEEASShellTangentDiagnostic : public FEDiagnostic
{
public:
    FEEASShellTangentDiagnostic(FEModel& fem);
    virtual ~FEEASShellTangentDiagnostic(){}
    
    FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
    bool Init();
    
    bool Run();
    
protected:
    void deriv_residual(matrix& ke);
    void print_matrix(matrix& m);
public:
    FEDiagnosticScenario* m_pscn;
};

#endif /* FEShellTangentDiagnostic_hpp */
