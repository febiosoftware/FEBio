//
//  FETiedBiphasicDiagnostic.hpp
//  FEBioTest
//
//  Created by Gerard Ateshian on 1/29/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FETiedBiphasicDiagnostic_hpp
#define FETiedBiphasicDiagnostic_hpp

#include "FEDiagnostic.h"
#include "NumCore/DenseMatrix.h"
#include "FECore/log.h"

//-----------------------------------------------------------------------------
class FETiedBiphasicScenario : public FEDiagnosticScenario
{
public:
    FETiedBiphasicScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }
    
public:
    double	m_dt;
};

//-----------------------------------------------------------------------------
class FETiedBiphasicTangentHex8 : public FETiedBiphasicScenario
{
public:
    FETiedBiphasicTangentHex8(FEDiagnostic* pdia);
    
    bool Init() override;
    
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
class FETiedBiphasicTangentHex20 : public FETiedBiphasicScenario
{
public:
    FETiedBiphasicTangentHex20(FEDiagnostic* pdia);
    
    bool Init() override;
    
    DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
class FETiedBiphasicDiagnostic : public FEDiagnostic
{
public:
    FETiedBiphasicDiagnostic(FEModel& fem);
    ~FETiedBiphasicDiagnostic();
    
    bool Run();
    
    bool Init();
    
    FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void print_matrix(matrix& m);
    void print_matrix(SparseMatrix& m);
    
    void deriv_residual(matrix& K);
    
public:
    FETiedBiphasicScenario*  m_pscn;
};

#endif /* FETiedBiphasicDiagnostic_hpp */
