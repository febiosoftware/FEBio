#pragma once
#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FEFluidFSIScenario : public FEDiagnosticScenario
{
public:
    FEFluidFSIScenario(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_dt = 1.0; }
    
public:
    double	m_dt;
};

//-----------------------------------------------------------------------------
class FEFluidFSITangentUniaxial : public FEFluidFSIScenario
{
public:
    FEFluidFSITangentUniaxial(FEDiagnostic* pdia);
    
    bool Init() override;
    
private:
    double		m_dilation;
    
    DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! The FEFluidFSITangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FEFluidFSITangentDiagnostic : public FEDiagnostic
{
public:
    FEFluidFSITangentDiagnostic(FEModel& fem);
    virtual ~FEFluidFSITangentDiagnostic(){}
    
    bool Init();
    
    bool Run();
    
    FEDiagnosticScenario* CreateScenario(const std::string& sname);
    
protected:
    void deriv_residual(matrix& ke);
    
    void print_matrix(matrix& m);
    
public:
    FEFluidFSIScenario*	m_pscn;
};
