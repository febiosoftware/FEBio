#pragma once
#include "FEElasticSolidDomain.h"
#include "FECore\tens3d.h"

//-----------------------------------------------------------------------------
//! This class implements a domain used in an elastic remodeling problem.
//! It differs from the FEElasticSolidDomain in that it adds a stiffness matrix
//! due to the deformation dependent density.
class FEElasticDomain2O : public FEElasticSolidDomain
{
public:
	//! constructor
	FEElasticDomain2O(FEMesh* pm, FEMaterial* pmat);
	
	void InitElements();
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	void UpdateElementStress(int iel, double dt);

public:
	// --- S T I F F N E S S ---
	//! calculates the solid element stiffness matrix
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);
	void ElementMaterialStiffness(FESolidElement &el, matrix &ke);

	void defhess(FESolidElement &el, tens3drs &G, int n);
};
