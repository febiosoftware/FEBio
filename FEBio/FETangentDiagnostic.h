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
	FETangentDiagnostic(FEM& fem);
	virtual ~FETangentDiagnostic();

	bool Init();

	bool Run();

	bool ParseSection(XMLTag& tag);

protected:
	static double residual(double d);

protected:
	void solve();

	void deriv_residual(matrix& ke);

	double	m_strain;
	double	m_stretch;

	static FETangentDiagnostic* m_pthis;
};

#endif // !defined(AFX_FETANGENTDIAGNOSTIC_H__C41CFF58_F916_4835_9993_7B461D9F282B__INCLUDED_)
