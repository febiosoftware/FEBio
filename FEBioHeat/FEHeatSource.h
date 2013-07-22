#pragma once
#include <FECore/FEParameterList.h>
#include <FECore/FEGlobalVector.h>
#include <FECore/FEBodyForce.h>
#include "FEHeatSolidDomain.h"

//-----------------------------------------------------------------------------
// Forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
class FEHeatSource : public FEBodyLoad
{
public:
	//! constructor
	FEHeatSource(FEModel* pfem);

	//! calculate the RHS-contribution of the source
	void Residual(FEGlobalVector& R);

protected:
	void ElementResidual(FEHeatSolidDomain& dom, FESolidElement& el, vector<double>& fe);

protected:
	FEModel*	m_pfem;		//!< the model this heat source belongs to

public:
	double	m_Q;	// source value

	DECLARE_PARAMETER_LIST();
};
