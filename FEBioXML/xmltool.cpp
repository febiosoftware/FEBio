#include "stdafx.h"
#include "xmltool.h"


//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool fexml::readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName)
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
        default: break;
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void fexml::readList(XMLTag& tag, vector<int>& l)
{
	// make sure the list is empty
	l.clear();

	// get a pointer to the value
	const char* sz = tag.szvalue();

	// parse the string
	const char* ch;
	do
	{
		int n0, n1, nn;
		int nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
			break;
		}

		for (int i = n0; i <= n1; i += nn) l.push_back(i);

		ch = strchr(sz, ',');
		if (ch) sz = ch + 1;
	}
	while (ch != 0);
}
