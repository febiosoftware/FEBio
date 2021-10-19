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
#include "FEBioImport.h"

class FEFacetSet;

//-----------------------------------------------------------------------------
// Constraints Section
// (Base class. Don't use this directly!)
class FEBioConstraintsSection : public FEFileSection
{
public:
	FEBioConstraintsSection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	bool ParseSurfaceSection(XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 1.x)
class FEBioConstraintsSection1x : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection1x(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidConstraint(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.0)
class FEBioConstraintsSection2 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection2(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
	
protected:
	void ParseRigidConstraint20(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Constraints Section (format 2.5)
class FEBioConstraintsSection25 : public FEBioConstraintsSection
{
public:
	FEBioConstraintsSection25(FEFileImport* pim) : FEBioConstraintsSection(pim){}
	void Parse(XMLTag& tag);
};
