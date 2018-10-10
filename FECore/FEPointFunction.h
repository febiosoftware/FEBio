#pragma once
#include "FEFunction1D.h"

#include <vector>

//-----------------------------------------------------------------------------
class DumpStream;

//-----------------------------------------------------------------------------
//! This class implements a function defined by a set of ordered (point,value) pairs.
//! It uses an interpolation scheme to interpolate between data points and also 
//! has a way of specifying the function value outside of the domain defined by the first and 
//! last data point. 

class FECORE_API FEPointFunction : public FEFunction1D
{
public:
	//! Load point structure
	struct LOADPOINT
	{
		double time;
		double value;
	};

public:
	//! Interpolation functions
	enum INTFUNC { STEP = 0, LINEAR = 1, SMOOTH = 2 };

	//! Extend mode
	enum EXTMODE { CONSTANT, EXTRAPOLATE, REPEAT, REPEAT_OFFSET };

public:
	//! default constructor
	FEPointFunction(FEModel* fem);

	//! destructor
	virtual ~FEPointFunction() {}

	//! adds a point to the point curve
	void Add(double x, double y);

	//! Clears the loadcurve data
	void Clear();

	//! set the x and y value of point i
	void SetPoint(int i, double x, double y);

	//! Set the type of interpolation
	void SetInterpolation(INTFUNC fnc) { m_fnc = fnc; }

	//! Set the extend mode
	void SetExtendMode(EXTMODE mode) { m_ext = mode; }

	//! returns point i
	LOADPOINT LoadPoint(int i) const;

	//! finds closest load point
	int FindPoint(double t, double& tval, int startIndex = 0);

	//! return nr of points
	int Points() const;

	//! see if there is a point at time t
	bool HasPoint(double t) const;

	//! Serialize data to archive
	void Serialize(DumpStream& ar);

	// copy data from other curve
	FEFunction1D* copy();

public: // implement from base class

		//! returns the value of the load curve at time
	double value(double x) const;

	//! returns the derivative value at time
	double derive(double x) const;

protected:
	double ExtendValue(double t) const;

protected:
	std::vector<vec2d>	m_points;

	int		m_fnc;	//!< interpolation function
	int		m_ext;	//!< extend mode

	DECLARE_FECORE_CLASS();
};

typedef FEPointFunction::LOADPOINT LOADPOINT;
