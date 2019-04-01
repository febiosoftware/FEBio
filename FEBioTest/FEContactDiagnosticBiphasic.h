#pragma once
#include "FEDiagnostic.h"
#include <FECore/DenseMatrix.h>

//-----------------------------------------------------------------------------
class FEContactBiphasicScenario : public FEDiagnosticScenario
{
public:
    FEContactBiphasicScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }
    
public:
    double	m_dt;
};

//-----------------------------------------------------------------------------
class FEContactBiphasicTangentHex8 : public FEContactBiphasicScenario
{
public:
    FEContactBiphasicTangentHex8(FEDiagnostic* pdia);
    
    bool Init() override;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEContactBiphasicTangentHex20 : public FEContactBiphasicScenario
{
public:
    FEContactBiphasicTangentHex20(FEDiagnostic* pdia);
    
    bool Init() override;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEContactDiagnosticBiphasic : public FEDiagnostic
{
public:
    FEContactDiagnosticBiphasic(FEModel& fem);
    ~FEContactDiagnosticBiphasic();
    
    bool Run();
    
    bool Init();
    
    FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void print_matrix(matrix& m);
    void print_matrix(SparseMatrix& m);
    
    void deriv_residual(matrix& K);
    
public:
    FEContactBiphasicScenario*  m_pscn;
};
