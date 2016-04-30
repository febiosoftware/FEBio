#pragma once

//-----------------------------------------------------------------------------
class FEParam;

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
class FEParamValidator
{
public:
	FEParamValidator(){}
	virtual ~FEParamValidator(){}

	//! Function that is used to see if a parameter is valid
	virtual bool is_valid(const FEParam& p) const = 0;

	//! Creates a copy of the validator
	virtual FEParamValidator* copy() const = 0;
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

private:
	FEParamRange	m_rng;
	int				m_nmin;
	int				m_nmax;
};

//-----------------------------------------------------------------------------
//! Validator for FE_PARAM_DOUBLE parameters that does basic range checking.
class FEDoubleValidator : public FEParamValidator_<FEDoubleValidator>
{
public:
	FEDoubleValidator(FEParamRange rng, double fmin, double fmax) : m_rng(rng), m_fmin(fmin), m_fmax(fmax) {}

	//! See if the parameter is an FE_PARAM_DOUBLE and within the specified range
	bool is_valid(const FEParam& p) const;

private:
	FEParamRange	m_rng;
	double			m_fmin;
	double			m_fmax;
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
