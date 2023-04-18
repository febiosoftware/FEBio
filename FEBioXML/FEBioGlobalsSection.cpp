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
#include "FEBioGlobalsSection.h"
#include "FECore/FEModel.h"
#include "FECore/FEGlobalData.h"
#include "FECore/FECoreKernel.h"

//-----------------------------------------------------------------------------
//!  This function reads the global variables from the xml file
//!
void FEBioGlobalsSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if      (tag == "Constants"          ) ParseConstants(tag);
		else if (tag == "Solutes"            ) ParseGlobalData(tag);
		else if (tag == "SolidBoundMolecules") ParseGlobalData(tag);
		else if (tag == "Variables"          ) ParseVariables(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseConstants(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	++tag;
	string s;
	double v;
	do
	{
		s = string(tag.Name());
		tag.value(v);
		fem.SetGlobalConstant(s, v);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseGlobalData(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	
	// read the global solute data
	++tag;
	do
	{
		// create new global data
		FEGlobalData* pgd = fecore_new<FEGlobalData>(tag.Name(), &fem);
		if (pgd == 0) throw XMLReader::InvalidTag(tag);
		fem.AddGlobalData(pgd);

		// TODO: We have to call the Init member here because solute data 
		//       allocates the concentration dofs and they have to be allocated before 
		//       materials are read in. I'd like to move this to FEModel::Init but not sure
		//       yet how. 
		pgd->Init();

		// read solute properties
		ReadParameterList(tag, pgd);
		
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioGlobalsSection::ParseVariables(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	fem.GetParameterList();
	if (tag.isleaf()) return;
	++tag;
	do
	{
		if (tag == "var")
		{
			const char* szname = tag.AttributeValue("name");
			double v = 0; tag.value(v);
			fem.AddGlobalVariable(szname, v);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}
