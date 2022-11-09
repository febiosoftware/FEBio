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
#include "FEBioModuleSection.h"

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//! It is currently only used to allocate the FESolver.
void FEBioModuleSection::Parse(XMLTag &tag)
{
	int nversion = GetFileReader()->GetFileVersion();

	// For version 2.5 and up this tag can only be read in once and we 
	// determine that by inspecting the m_szmod member.
	if (nversion >= 0x0205)
	{
		std::string moduleName = GetBuilder()->GetModuleName();
		if (moduleName.empty() == false) throw XMLReader::InvalidTag(tag);
	}

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	// some special case
	if (strcmp(szt, "explicit-solid") == 0)
	{
		szt = "solid";
		GetBuilder()->SetDefaultSolver("explicit-solid");
	}
	else if (strcmp(szt, "CG-solid") == 0)
	{
		szt = "solid";
		GetBuilder()->SetDefaultSolver("CG-solid");
	}

	GetBuilder()->SetActiveModule(szt);

	if (tag.isempty()) return;

	++tag;
	do {
		if (tag == "units")
		{
			const char* szunits = tag.szvalue();
			GetFEModel()->SetUnits(szunits);
		}
		
		++tag;
	} while (!tag.isend());
}
