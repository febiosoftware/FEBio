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
//! The FETangentDiagnostic class tests the stiffness matrix implementation
//! by comparing it to a numerical approximating of the derivative of the
//! residual.

class FETangentDiagnostic : public FEDiagnostic  
{
public:
	enum TD_Scenario {
		TDS_UNIAXIAL,
		TDS_SIMPLE_SHEAR
	};

public:
	FETangentDiagnostic(FEM& fem);
	virtual ~FETangentDiagnostic(){}

	bool Init();

	bool Run();

protected:
	void BuildUniaxial();
	void BuildSimpleShear();

	void deriv_residual(matrix& ke);

public:
	TD_Scenario	m_scn;
	double		m_strain;
};

#endif // !defined(AFX_FETANGENTDIAGNOSTIC_H__C41CFF58_F916_4835_9993_7B461D9F282B__INCLUDED_)
