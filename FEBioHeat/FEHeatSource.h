#pragma once
#include "FECore/FEParameterList.h"
#include "FECore/FEGlobalVector.h"
#include "FECore/FEBodyLoad.h"

//-----------------------------------------------------------------------------
// Forward declaration of FEModel class
class FEModel;

//-----------------------------------------------------------------------------
// Base class for heat sources
class FEHeatSource : public FEBodyLoad
{
public:
	//! constructor
	FEHeatSource(FEModel* pfem);

public:
	virtual double value(FEMaterialPoint& mp) = 0;
};

//-----------------------------------------------------------------------------
class FEConstHeatSource : public FEHeatSource
{
public:
	FEConstHeatSource(FEModel* fem);

	double value(FEMaterialPoint& mp) override { return m_Q; }

private:
	double	m_Q;	// source value

	DECLARE_PARAMETER_LIST();
};
