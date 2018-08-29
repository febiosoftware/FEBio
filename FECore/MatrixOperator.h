#pragma once
#include "fecore_api.h"

// abstract base class for matrix operators, i.e. a class that can calculate a matrix-vector product
class FECORE_API MatrixOperator
{
public:
	MatrixOperator() {}
	virtual ~MatrixOperator() {}

	// calculate the product Ax = y
	virtual bool mult_vector(double* x, double* y) = 0;
};
