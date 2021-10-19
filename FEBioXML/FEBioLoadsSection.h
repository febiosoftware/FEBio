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

//-----------------------------------------------------------------------------
// Version 1.2
class FEBioLoadsSection1x : public FEFileSection
{
public:
	FEBioLoadsSection1x(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseBodyForce  (XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Version 2.0
class FEBioLoadsSection2 : public FEFileSection
{
public:
	FEBioLoadsSection2(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
	void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
	void ParseSurfaceLoadSurface(XMLTag& tag, FESurface* psurf);
};

//-----------------------------------------------------------------------------
// Version 2.5
class FEBioLoadsSection25 : public FEFileSection
{
public:
	FEBioLoadsSection25(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);

protected:
	void ParseObsoleteLoad(XMLTag& tag);
};

//-----------------------------------------------------------------------------
// Version 3.0
class FEBioLoadsSection3 : public FEFileSection
{
public:
	FEBioLoadsSection3(FEFileImport* pim) : FEFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalLoad  (XMLTag& tag);
	void ParseEdgeLoad   (XMLTag& tag);
	void ParseSurfaceLoad(XMLTag& tag);
	void ParseBodyLoad   (XMLTag& tag);
};
