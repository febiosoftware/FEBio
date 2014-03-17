#pragma once
#include "FECore/FEParameterList.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/FEBodyLoad.h"
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

	//! serialization
	void Serialize(DumpFile& ar);

protected:
	void ElementResidual(FEHeatSolidDomain& dom, FESolidElement& el, vector<double>& fe);

public:
	double	m_Q;	// source value

	DECLARE_PARAMETER_LIST();
};
