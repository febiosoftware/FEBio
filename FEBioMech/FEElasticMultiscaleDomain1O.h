#pragma once
#include "FEElasticSolidDomain.h"
#include "FECore/tens3d.h"

//-----------------------------------------------------------------------------
//! This class implements a domain used in an elastic remodeling problem.
//! It differs from the FEElasticSolidDomain in that it adds a stiffness matrix
//! due to the deformation dependent density.
class FEElasticMultiscaleDomain1O : public FEElasticSolidDomain
{
public:
	//! constructor
	FEElasticMultiscaleDomain1O(FEModel* pfem);
	
	void InitElements();

	void UpdateElementStress(int iel, double dt);

};
