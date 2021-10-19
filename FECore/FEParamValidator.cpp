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



#include "stdafx.h"
#include "FEParamValidator.h"
#include "FEParam.h"
#include "FECoreKernel.h"
#include "DumpStream.h"
#include "FEModelParam.h"

//-----------------------------------------------------------------------------
bool is_inside_range_int(int ival, int rng, int imin, int imax)
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
bool is_inside_range_double(double val, int rng, double dmin, double dmax)
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

	return bvalid;
}

//-----------------------------------------------------------------------------
void FEIntValidator::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_rng;
	ar & m_nmin & m_nmax;
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

	return bvalid;
}

//-----------------------------------------------------------------------------
void FEDoubleValidator::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_rng;
	ar & m_fmin & m_fmax;
}

//-----------------------------------------------------------------------------
bool FEParamDoubleValidator::is_valid(const FEParam& p) const
{
	if (p.type() != FE_PARAM_DOUBLE_MAPPED) return false;
	const FEParamDouble& d = p.value<FEParamDouble>();

	// This only works for const values
	double val = 0.0;
	bool bvalid = true;
	if (d.isConst())
	{
		val = d.constValue();
		bvalid = is_inside_range_double(val, m_rng, m_fmin, m_fmax);
	}

	return bvalid;
}

//-----------------------------------------------------------------------------
void FEParamDoubleValidator::Serialize(DumpStream& ar)
{
	if (ar.IsShallow()) return;
	ar & m_rng;
	ar & m_fmin & m_fmax;
}
