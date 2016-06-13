#pragma once
#include "FEElasticSolidDomain.h"

//-----------------------------------------------------------------------------
//! This class implements a domain used in an elastic remodeling problem.
//! It differs from the FEElasticSolidDomain in that it adds a stiffness matrix
//! due to the deformation dependent density.
class FERemodelingElasticDomain : public FEElasticSolidDomain
{
public:
	//! constructor
	FERemodelingElasticDomain(FEModel* pfem);

	//! reset element data
	void Reset();

	//! initialize class
	bool Initialize();

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver);

	//! calculates the solid element stiffness matrix (\todo is this actually used anywhere?)
	virtual void ElementStiffness(const FETimeInfo& tp, int iel, matrix& ke);

private:
	//! density stiffness component
	void ElementDensityStiffness(double dt, FESolidElement& el, matrix& ke);
};
