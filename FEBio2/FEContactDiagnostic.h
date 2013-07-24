// FEContactDiagnostic.h: interface for the FEContactDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FECONTACTDIAGNOSTIC_H__6DAAF8C6_9D08_49F7_955B_08201F1A48C9__INCLUDED_)
#define AFX_FECONTACTDIAGNOSTIC_H__6DAAF8C6_9D08_49F7_955B_08201F1A48C9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FEDiagnostic.h"
#include "NumCore/DenseMatrix.h"

class FEContactDiagnostic : public FEDiagnostic  
{
public:
	FEContactDiagnostic(FEModel& fem);
	virtual ~FEContactDiagnostic();

	bool Run();

	bool Init();

protected:
	void deriv_residual(DenseMatrix& K);
};

#endif // !defined(AFX_FECONTACTDIAGNOSTIC_H__6DAAF8C6_9D08_49F7_955B_08201F1A48C9__INCLUDED_)
