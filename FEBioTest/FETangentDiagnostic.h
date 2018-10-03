// FETangentDiagnostic.h: interface for the FETangentDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FETANGENTDIAGNOSTIC_H__C41CFF58_F916_4835_9993_7B461D9F282B__INCLUDED_)
#define AFX_FETANGENTDIAGNOSTIC_H__C41CFF58_F916_4835_9993_7B461D9F282B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

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

public:
	FEDiagnosticScenario* m_pscn;
};

#endif // !defined(AFX_FETANGENTDIAGNOSTIC_H__C41CFF58_F916_4835_9993_7B461D9F282B__INCLUDED_)
