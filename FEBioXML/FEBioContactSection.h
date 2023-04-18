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
#include "FileImport.h"
#include <FECore/FESurfacePairConstraint.h>

//-----------------------------------------------------------------------------
// Contact section (new in version 2.0)
class FEBioContactSection : public FEFileSection
{
protected:
	//! missing primary surface
	class MissingPrimarySurface : public FEFileException
	{
	public:
		MissingPrimarySurface();
	};

	//! missing secondary surface
	class MissingSecondarySurface : public FEFileException
	{
	public:
		MissingSecondarySurface();
	};

public:
	FEBioContactSection(FEFileImport* pim) : FEFileSection(pim){}

protected:
	void ParseLinearConstraint     (XMLTag& tag);

protected:
	bool ParseSurfaceSection  (XMLTag& tag, FESurface& s, int nfmt, bool bnodal);
};

//-----------------------------------------------------------------------------
// Version 2.0
class FEBioContactSection2 : public FEBioContactSection
{
public:
	FEBioContactSection2(FEFileImport* im) : FEBioContactSection(im){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidInterface(XMLTag& tag);
	void ParseRigidWall(XMLTag& tag);
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci);
};

//-----------------------------------------------------------------------------
// Version 2.5
class FEBioContactSection25 : public FEBioContactSection
{
public:
	FEBioContactSection25(FEFileImport* im) : FEBioContactSection(im){}
	void Parse(XMLTag& tag);

protected:
	void ParseRigidWall(XMLTag& tag);
	void ParseRigidSliding(XMLTag& tag);
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci);
};


//-----------------------------------------------------------------------------
// Version 4.0
class FEBioContactSection4 : public FEBioContactSection
{
public:
	FEBioContactSection4(FEFileImport* im) : FEBioContactSection(im) {}
	void Parse(XMLTag& tag);

protected:
	void ParseContactInterface(XMLTag& tag, FESurfacePairConstraint* pci);
};
