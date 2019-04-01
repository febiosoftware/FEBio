#pragma once
#include "FEDiagnostic.h"
#include <FECore/DenseMatrix.h>

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
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FETiedBiphasicTangentHex20 : public FETiedBiphasicScenario
{
public:
    FETiedBiphasicTangentHex20(FEDiagnostic* pdia);
    
    bool Init() override;
    
    DECLARE_FECORE_CLASS();
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
