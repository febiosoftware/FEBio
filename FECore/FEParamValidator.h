/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class FEParam;
class DumpStream;

//-----------------------------------------------------------------------------
//! Range for parameters
//! Used only by parameters of type FE_PARAM_INT and FE_PARAM_DOUBLE
enum FEParamRange {
	FE_DONT_CARE,			// parameter can take on any value
	FE_GREATER,				// parameter must be greater than val (> val)
	FE_GREATER_OR_EQUAL,	// parameter must be greater or equal than val ( >= val)
	FE_LESS,				// parameter must be less than val ( < val)
	FE_LESS_OR_EQUAL,		// parameter must be less or equal than val ( <= val)
	FE_OPEN,				// parameter must be inside the open interval (min,max)
	FE_CLOSED,				// parameter must be inside the closed interval [min, max]
	FE_LEFT_OPEN,			// parameter must be inside the half-open interval (min, max]
	FE_RIGHT_OPEN,			// parameter must be inside the right-open interval [min, max)
	FE_NOT_EQUAL,			// parameter must not equal val ( != val)
};

//-----------------------------------------------------------------------------
//! Abstract base class for parameter validators. 
class FECORE_API FEParamValidator
{
public:
	FEParamValidator(){}
	virtual ~FEParamValidator(){}

	//! Function that is used to see if a parameter is valid
	virtual bool is_valid(const FEParam& p) const = 0;

	//! Creates a copy of the validator
	virtual FEParamValidator* copy() const = 0;

	//! Serialize validator data
	virtual void Serialize(DumpStream& ar) = 0;
};

//-----------------------------------------------------------------------------
//! Base class for parameter validators that uses CRTP for implementing a copy method.
//! Deriving validators from this class automatically adds the copy method which simplifies
//! the implementation of new validators. It does assume that the derived class has a valid
//! copy constructor.
template <class T> class FEParamValidator_ : public FEParamValidator
{
public:
	FEParamValidator* copy() const { return new T(static_cast<T const&>(*this)); }
};

//-----------------------------------------------------------------------------
//! Validator for FE_PARAM_INT parameters that does basic range checking.
class FEIntValidator : public FEParamValidator_<FEIntValidator>
{
public:
	FEIntValidator(FEParamRange rng, int nmin, int nmax) : m_rng(rng), m_nmin(nmin), m_nmax(nmax) {}

	//! See if the parameter is an FE_PARAM_INT and within the specified range
	bool is_valid(const FEParam& p) const;

	void Serialize(DumpStream& ar);

private:
	int		m_rng;
	int		m_nmin;
	int		m_nmax;
};

//-----------------------------------------------------------------------------
//! Validator for FE_PARAM_DOUBLE parameters that does basic range checking.
class FEDoubleValidator : public FEParamValidator_<FEDoubleValidator>
{
public:
	FEDoubleValidator(FEParamRange rng, double fmin, double fmax) : m_rng(rng), m_fmin(fmin), m_fmax(fmax) {}

	//! See if the parameter is an FE_PARAM_DOUBLE and within the specified range
	bool is_valid(const FEParam& p) const;

	void Serialize(DumpStream& ar);

private:
	int			m_rng;
	double		m_fmin;
	double		m_fmax;
};

//-----------------------------------------------------------------------------
//! Validator for FE_PARAM_DOUBLE parameters that does basic range checking.
class FEParamDoubleValidator : public FEParamValidator_<FEParamDoubleValidator>
{
public:
	FEParamDoubleValidator(FEParamRange rng, double fmin, double fmax) : m_rng(rng), m_fmin(fmin), m_fmax(fmax) {}

	//! See if the parameter is an FE_PARAM_DOUBLE and within the specified range
	bool is_valid(const FEParam& p) const;

	void Serialize(DumpStream& ar);

private:
	int			m_rng;
	double		m_fmin;
	double		m_fmax;
};


//-----------------------------------------------------------------------------
// helper class for defining ranges
class RANGE
{
public:
	FEParamRange	m_rt;
	double			m_fmin;
	double			m_fmax;
public:
	RANGE(FEParamRange rt) { m_rt = rt; m_fmin = m_fmax = 0.0; }
	RANGE(FEParamRange rt, double fval) { m_rt = rt; m_fmin = fval; m_fmax = 0.0; }
	RANGE(FEParamRange rt, double fmin, double fmax) { m_rt = rt; m_fmin = fmin; m_fmax = fmax; }
};

inline RANGE FE_RANGE_DONT_CARE() { return RANGE(FE_DONT_CARE); }
inline RANGE FE_RANGE_GREATER         (double f) { return RANGE(FE_GREATER         , f); }
inline RANGE FE_RANGE_GREATER_OR_EQUAL(double f) { return RANGE(FE_GREATER_OR_EQUAL, f); }
inline RANGE FE_RANGE_LESS            (double f) { return RANGE(FE_LESS            , f); }
inline RANGE FE_RANGE_LESS_OR_EQUAL   (double f) { return RANGE(FE_LESS_OR_EQUAL   , f); }
inline RANGE FE_RANGE_OPEN      (double fmin, double fmax) {return RANGE(FE_OPEN      , fmin, fmax); }
inline RANGE FE_RANGE_CLOSED    (double fmin, double fmax) {return RANGE(FE_CLOSED    , fmin, fmax); }
inline RANGE FE_RANGE_LEFT_OPEN (double fmin, double fmax) {return RANGE(FE_LEFT_OPEN , fmin, fmax); }
inline RANGE FE_RANGE_RIGHT_OPEN(double fmin, double fmax) {return RANGE(FE_RIGHT_OPEN, fmin, fmax); }
inline RANGE FE_RANGE_NOT_EQUAL (double f) { return RANGE(FE_NOT_EQUAL, f); }
