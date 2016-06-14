#pragma once
#include "FEElasticSolidDomain2O.h"
#include <FECore/tens3d.h>
#include <FECore/tens5d.h>
#include <FECore/tens6d.h>
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
//! This class implements a domain used in an elastic remodeling problem.
//! It differs from the FEElasticSolidDomain in that it adds a stiffness matrix
//! due to the deformation dependent density.
class FEElasticMultiscaleDomain2O : public FEElasticSolidDomain2O
{
public:
	//! constructor
	FEElasticMultiscaleDomain2O(FEModel* pfem);
	
	//! initialize class
	bool Initialize();

	//! Update 
	void Update(const FETimeInfo& timeInfo);
};
