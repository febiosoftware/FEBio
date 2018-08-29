#pragma once
#include "LoadCurve.h"

#include <vector>

//-----------------------------------------------------------------------------
class DumpStream;

//-----------------------------------------------------------------------------
//! This class implements the concept of a loadcurve.

//! A loadcurve is basically a discretized function of time versus load,
//! where load can be interpreted differently in different contexts.
//! The loadcurve stores the (time,load)-pairs for fixed points, which
//! are input from the input file.
//! In between timesteps, the loadcurve class interpolates the (time,load)
//! data pairs according to the interpolation function.

class FECORE_API FEDataLoadCurve : public FELoadCurve
{
public:
	class FEDataPoint : public FECoreBase
	{
	public:
		FEDataPoint() : FECoreBase(FEOBJECT_ID), x(0), y(0) {}
		FEDataPoint(double X, double Y) : FECoreBase(FEOBJECT_ID), x(X), y(Y) {}

	public:
		double	x, y;

		DECLARE_PARAMETER_LIST();
	};

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
	FEDataLoadCurve(FEModel* fem);

	//! destructor
	virtual ~FEDataLoadCurve() {}

	//! adds a point to the loadcurve
	void Add(double time, double value);

	//! Clears the loadcurve data
	void Clear();

	//! set the time and data value of point i of the load curve
	void SetPoint(int i, double time, double val);

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
	bool CopyFrom(FELoadCurve* lc);

public: // implement from base class

	//! returns the value of the load curve at time
	double Value(double time) const;

	//! returns the derivative value at time
	double Deriv(double time) const;

protected:
	double ExtendValue(double t) const;

protected:
	FEVecPropertyT<FEDataPoint>	m_points;

	INTFUNC		m_fnc;	//!< interpolation function
	EXTMODE		m_ext;	//!< extend mode
};

typedef FEDataLoadCurve::LOADPOINT LOADPOINT;
