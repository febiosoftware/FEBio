#include "stdafx.h"
#include "FEParamValidator.h"
#include "FEParam.h"
#include "FECoreKernel.h"

//-----------------------------------------------------------------------------
bool is_inside_range_int(int ival, FEParamRange rng, int imin, int imax)
{
	switch (rng)
	{
	case FE_GREATER         : return (ival >  imin); break;
	case FE_GREATER_OR_EQUAL: return (ival >= imin); break;
	case FE_LESS            : return (ival <  imin); break;
	case FE_LESS_OR_EQUAL   : return (ival <= imin); break;
	case FE_OPEN            : return ((ival >  imin) && (ival <  imax)); break;
	case FE_CLOSED          : return ((ival >= imin) && (ival <= imax)); break;
	case FE_LEFT_OPEN       : return ((ival >  imin) && (ival <= imax)); break;
	case FE_RIGHT_OPEN      : return ((ival >= imin) && (ival <  imax)); break;
	case FE_NOT_EQUAL       : return (ival != imin); break;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool is_inside_range_double(double val, FEParamRange rng, double dmin, double dmax)
{
	switch (rng)
	{
	case FE_GREATER         : return (val >  dmin); break;
	case FE_GREATER_OR_EQUAL: return (val >= dmin); break;
	case FE_LESS            : return (val <  dmin); break;
	case FE_LESS_OR_EQUAL   : return (val <= dmin); break;
	case FE_OPEN            : return ((val >  dmin) && (val <  dmax)); break;
	case FE_CLOSED          : return ((val >= dmin) && (val <= dmax)); break;
	case FE_LEFT_OPEN       : return ((val >  dmin) && (val <= dmax)); break;
	case FE_RIGHT_OPEN      : return ((val >= dmin) && (val <  dmax)); break;
	case FE_NOT_EQUAL       : return (val != dmin); break;
	}
	return false;
}

//-----------------------------------------------------------------------------
bool FEIntValidator::is_valid(const FEParam& p) const
{
	if (p.type() != FE_PARAM_INT) return false;

	bool bvalid = true;
	int val = 0;
	if (p.dim() == 1)
	{
		val = p.value<int>();
		bvalid = is_inside_range_int(val, m_rng, m_nmin, m_nmax);
	}
	else 
	{
		for (int i = 0; i<p.dim(); ++i)
		{
			val = p.value<int>(i);
			bvalid = is_inside_range_int(val, m_rng, m_nmin, m_nmax);
			if (bvalid == false) break;
		}
	}

	if (bvalid == false)
	{
		char szerr[256] = {0};
		switch (m_rng)
		{
		case FE_GREATER         : sprintf(szerr, "%s (=%d) must be greater than %d"                    , p.name(), val, m_nmin); break;
		case FE_GREATER_OR_EQUAL: sprintf(szerr, "%s (=%d) must be greater than or equal to %d"        , p.name(), val, m_nmin); break;
		case FE_LESS            : sprintf(szerr, "%s (=%d) must be less than %d"                       , p.name(), val, m_nmin); break;
		case FE_LESS_OR_EQUAL   : sprintf(szerr, "%s (=%d) must be less than or equal to %d"           , p.name(), val, m_nmin); break;
		case FE_OPEN            : sprintf(szerr, "%s (=%d) must be in the open interval (%d, %d)"      , p.name(), val, m_nmin, m_nmax); break;
		case FE_CLOSED          : sprintf(szerr, "%s (=%d) must be in the closed interval [%d, %d]"    , p.name(), val, m_nmin, m_nmax); break;
		case FE_LEFT_OPEN       : sprintf(szerr, "%s (=%d) must be in the left-open interval (%d, %d]" , p.name(), val, m_nmin, m_nmax); break;
		case FE_RIGHT_OPEN      : sprintf(szerr, "%s (=%d) must be in the right-open interval [%d, %d)", p.name(), val, m_nmin, m_nmax); break;
		case FE_NOT_EQUAL       : sprintf(szerr, "%s (=%d) must not equal %d"                          , p.name(), m_nmin);
		default:
			sprintf(szerr, "%s has an invalid range");
		}
		
		return fecore_error(szerr);
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEDoubleValidator::is_valid(const FEParam& p) const
{
	if (p.type() != FE_PARAM_DOUBLE) return false;

	bool bvalid;
	double val = 0;
	if (p.dim() == 1)
	{
		val = p.value<double>();
		bvalid = is_inside_range_double(val, m_rng, m_fmin, m_fmax);
	}
	else
	{
		for (int i = 0; i<p.dim(); ++i)
		{
			val = p.value<double>(i);
			bvalid = is_inside_range_double(val, m_rng, m_fmin, m_fmax);
			if (bvalid == false) break;
		}
	}

	if (bvalid == false)
	{
		char szerr[256] = { 0 };
		switch (m_rng)
		{
		case FE_GREATER         : sprintf(szerr, "%s (=%lg) must be greater than %lg"                     , p.name(), val, m_fmin); break;
		case FE_GREATER_OR_EQUAL: sprintf(szerr, "%s (=%lg) must be greater than or equal to %lg"         , p.name(), val, m_fmin); break;
		case FE_LESS            : sprintf(szerr, "%s (=%lg) must be less than %lg"                        , p.name(), val, m_fmin); break;
		case FE_LESS_OR_EQUAL   : sprintf(szerr, "%s (=%lg) must be less than or equal to %lg"            , p.name(), val, m_fmin); break;
		case FE_OPEN            : sprintf(szerr, "%s (=%lg) must be in the open interval (%lg, %lg)"      , p.name(), val, m_fmin, m_fmax); break;
		case FE_CLOSED          : sprintf(szerr, "%s (=%lg) must be in the closed interval [%lg, %lg]"    , p.name(), val, m_fmin, m_fmax); break;
		case FE_LEFT_OPEN       : sprintf(szerr, "%s (=%lg) must be in the left-open interval (%lg, %lg]" , p.name(), val, m_fmin, m_fmax); break;
		case FE_RIGHT_OPEN      : sprintf(szerr, "%s (=%lg) must be in the right-open interval [%lg, %lg)", p.name(), val, m_fmin, m_fmax); break;
		case FE_NOT_EQUAL       : sprintf(szerr, "%s (=%lg) must not equal %lg"                           , p.name(), m_fmin);
		default:
			sprintf(szerr, "%s has an invalid range");
		}
		return fecore_error(szerr);
	}

	return true;
}
