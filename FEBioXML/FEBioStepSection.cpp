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
#include "FEBioStepSection.h"
#include "FEBioModuleSection.h"
#include "FEBioControlSection.h"
#include "FEBioConstraintsSection.h"
#include "FEBioBoundarySection.h"
#include "FEBioLoadsSection.h"
#include "FEBioContactSection.h"
#include "FEBioInitialSection.h"
#include "FEBioBoundarySection3.h"

//-----------------------------------------------------------------------------
void FEBioStepSection::Parse(XMLTag& tag)
{
	// create next step
	FEBioImport* feb = GetFEBioImport();
	feb->GetBuilder()->NextStep();

	FEFileSectionMap Map;
	Map["Module"     ] = new FEBioModuleSection       (feb);
	Map["Control"    ] = new FEBioControlSection      (feb);
	Map["Constraints"] = new FEBioConstraintsSection1x(feb);
	Map["Boundary"   ] = new FEBioBoundarySection1x   (feb);
	Map["Loads"      ] = new FEBioLoadsSection1x      (feb);
	Map["Initial"    ] = new FEBioInitialSection      (feb);

	// parse the file sections
	Map.Parse(tag);
}

//-----------------------------------------------------------------------------
void FEBioStepSection2::Parse(XMLTag& tag)
{
	// create next step
	FEBioImport* feb = GetFEBioImport();
	feb->GetBuilder()->NextStep();

	FEFileSectionMap Map;
	Map["Module"     ] = new FEBioModuleSection      (feb);
	Map["Control"    ] = new FEBioControlSection     (feb);
	Map["Constraints"] = new FEBioConstraintsSection2(feb);
	Map["Boundary"   ] = new FEBioBoundarySection2   (feb);
	Map["Loads"      ] = new FEBioLoadsSection2      (feb);
	Map["Initial"    ] = new FEBioInitialSection     (feb);
	Map["Contact"    ] = new FEBioContactSection2    (feb);

	// parse the file sections
	Map.Parse(tag);
}

//-----------------------------------------------------------------------------
void FEBioStepSection25::Parse(XMLTag& tag)
{
	// create next step
	GetBuilder()->NextStep();

	FEFileImport* imp = GetFileReader();
	FEFileSectionMap Map;
	Map["Control"    ] = new FEStepControlSection     (imp);
	Map["Constraints"] = new FEBioConstraintsSection25(imp);
	Map["Boundary"   ] = new FEBioBoundarySection25   (imp);
	Map["Loads"      ] = new FEBioLoadsSection25      (imp);
	Map["Initial"    ] = new FEBioInitialSection25    (imp);
	Map["Contact"    ] = new FEBioContactSection25    (imp);

	// parse the file sections
	Map.Parse(tag);
}
