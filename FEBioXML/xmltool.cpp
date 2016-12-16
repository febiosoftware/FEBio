#include "stdafx.h"
#include "xmltool.h"


//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool xmlReadParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName)
{
	// see if we can find this parameter
	FEParam* pp = paramList.Find((paramName == 0 ? tag.Name() : paramName));
	if (pp == 0) return false;
	
	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE  : { double d; tag.value(d); pp->value<double  >() = d; } break;
		case FE_PARAM_INT     : { int n; tag.value(n); pp->value<int>() = n; } break;
		case FE_PARAM_BOOL    : { bool b; tag.value(b); pp->value<bool>() = b; } break;
		default:
			assert(false);
			return false;
		}
	}
	else
	{
		switch (pp->type())
		{
		case FE_PARAM_INT   : { vector<int> d(pp->dim()); tag.value(&d[0], pp->dim()); } break;
		case FE_PARAM_DOUBLE: { vector<double> d(pp->dim()); tag.value(&d[0], pp->dim()); } break;
		}
	}

	return true;
}
