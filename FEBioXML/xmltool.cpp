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
#include "xmltool.h"
#include <FECore/FECoreKernel.h>

int enumValue(const char* val, const char* szenum);
bool is_number(const char* sz);

//-----------------------------------------------------------------------------
bool parseEnumParam(FEParam* pp, const char* val)
{
	// get the enums
	const char* ch = pp->enums();
	if (ch == nullptr) return false;

	// special cases
	if (strncmp(ch, "@factory_list", 13) == 0)
	{
		int classID = atoi(ch + 14);

		FECoreKernel& fecore = FECoreKernel::GetInstance();
		for (int i = 0, n = 0; i < fecore.FactoryClasses(); ++i)
		{
			const FECoreFactory* fac = fecore.GetFactoryClass(i);
			if (fac->GetSuperClassID() == classID)
			{
				if (strcmp(fac->GetTypeStr(), val) == 0)
				{
					pp->value<int>() = n;
					return true;
				}
				n++;
			}
		}
		return false;
	}

	int n = enumValue(val, ch);
	if (n != -1) pp->value<int>() = n;
	else
	{
		// see if the value is an actual number
		if (is_number(val))
		{
			n = atoi(val);
			pp->value<int>() = n;
		}
	}

	return (n != -1);
}

//-----------------------------------------------------------------------------
//! This function parses a parameter list
bool fexml::readParameter(XMLTag& tag, FEParameterList& paramList, const char* paramName)
{
	// see if we can find this parameter
	FEParam* pp = paramList.FindFromName((paramName == 0 ? tag.Name() : paramName));
	if (pp == 0) return false;
	
	if (pp->dim() == 1)
	{
		switch (pp->type())
		{
		case FE_PARAM_DOUBLE  : { double d; tag.value(d); pp->value<double  >() = d; } break;
		case FE_PARAM_INT     : 
		{ 
			const char* szenum = pp->enums();
			if (szenum == 0)
			{
				int n;
				tag.value(n); pp->value<int>() = n;
			}
			else
			{
				bool bfound = parseEnumParam(pp, tag.szvalue());
				if (bfound == false) throw XMLReader::InvalidValue(tag);
			}
		}
		break;
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

bool fexml::readParameterList(XMLTag& tag, FECoreBase* pc)
{
	// make sure this tag has children
	if (tag.isleaf()) return true;

	// process the parameter lists
	++tag;
	do
	{
		if (readParameter(tag, pc) == false) return false;
		++tag;
	} 
	while (!tag.isend());

	return true;
}

bool fexml::readParameter(XMLTag& tag, FECoreBase* pc)
{
	FEParameterList& PL = pc->GetParameterList();
	if (readParameter(tag, PL) == false)
	{
		// see if this is a property
		// if we get here, the parameter is not found.
		// See if the parameter container has defined a property of this name
		int n = pc->FindPropertyIndex(tag.Name());
		if (n >= 0)
		{
			FEProperty* prop = pc->PropertyClass(n);
			const char* sztype = tag.AttributeValue("type");

			// try to allocate the class
			FECoreBase* pp = fecore_new<FECoreBase>(prop->GetSuperClassID(), sztype, pc->GetFEModel());
			if (pp == nullptr) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

			prop->SetProperty(pp);

			// read the property data
			readParameterList(tag, pp);
		}
		else throw XMLReader::InvalidTag(tag);
	}

	return true;
}

void readClassVariable(XMLTag& tag, FEClassDescriptor::ClassVariable* vars)
{
	if (tag.isleaf()) return;

	++tag;
	do {
		if (tag.isleaf())
		{
			const char* szname = tag.Name();

			// see if the type attribute is defined
			const char* sztype = tag.AttributeValue("type", true);
			if (sztype)
			{
				FEClassDescriptor::ClassVariable* child = new FEClassDescriptor::ClassVariable(szname, sztype);
				vars->AddVariable(child);
			}
			else
			{
				const char* szval = tag.szvalue();
				FEClassDescriptor::SimpleVariable* var = new FEClassDescriptor::SimpleVariable(szname, szval);
				vars->AddVariable(var);
			}
		}
		else
		{
			const char* szname = tag.Name();
			const char* sztype = tag.AttributeValue("type");

			FEClassDescriptor::ClassVariable* child = new FEClassDescriptor::ClassVariable(szname, sztype);
			vars->AddVariable(child);
			readClassVariable(tag, child);
		}
		++tag;
	}
	while (!tag.isend());
}

// create a class descriptor from the current tag
FEClassDescriptor* fexml::readParameterList(XMLTag& tag)
{
	const char* sztype = tag.AttributeValue("type");
	FEClassDescriptor* cd = new FEClassDescriptor(sztype);

	FEClassDescriptor::ClassVariable* root = cd->Root();
	readClassVariable(tag, root);

	return cd;
}
