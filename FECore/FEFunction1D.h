#pragma once

//-----------------------------------------------------------------------------
class FEModel;
class DumpStream;

//-----------------------------------------------------------------------------
// Class that represents a 1D function
// Currently, only function defined via load curves are defined. That's why
// we need to FEModel class. 
class FEFunction1D
{
public:
	FEFunction1D(FEModel* pfem);

	// Set the load curve index and scale factor
	void SetLoadCurveIndex(int lc, double scale = 1.0);

	// serialization
	void Serialize(DumpStream& ar);

public: 
	// evaluate the function at x
	double value(double x) const;

	// value of first derivative of function at x
	double derive(double x) const;

private:
	int		m_nlc;		//!< load curve index
	double	m_scale;	//!< load curve scale factor
	FEModel&	m_fem;
};
