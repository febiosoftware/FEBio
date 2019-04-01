#pragma once
#include "FEDiagnostic.h"

//-----------------------------------------------------------------------------
class FETangentUniaxial : public FEDiagnosticScenario
{
public:
	FETangentUniaxial(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_strain = 0.0; }

	bool Init() override;

private:
	double	m_strain;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FETangentSimpleShear : public FEDiagnosticScenario
{
public:
	FETangentSimpleShear(FEDiagnostic* pdia) : FEDiagnosticScenario(pdia) { m_strain = 0.0; }

	bool Init() override;

private:
	double	m_strain;

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! The FETangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FETangentDiagnostic : public FEDiagnostic  
{
public:
	FETangentDiagnostic(FEModel& fem);
	virtual ~FETangentDiagnostic(){}

	FEDiagnosticScenario* CreateScenario(const std::string& sname);

	bool Init();

	bool Run();

protected:
	void deriv_residual(matrix& ke);

	void print_matrix(matrix& m);

public:
	FEDiagnosticScenario* m_pscn;
};
