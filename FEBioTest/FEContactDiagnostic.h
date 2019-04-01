#pragma once
#include "FEDiagnostic.h"
#include <FECore/DenseMatrix.h>

class FEContactDiagnostic : public FEDiagnostic  
{
public:
	FEContactDiagnostic(FEModel& fem);
	virtual ~FEContactDiagnostic();

	bool Run();

	bool Init();

	void print_matrix(DenseMatrix& m);

protected:
	void deriv_residual(DenseMatrix& K);
};
