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
#include "FEBModel.h"

//-----------------------------------------------------------------------------
// Geometry Section (base class)
class FEBioGeometrySection : public FEBioFileSection
{
public:
	FEBioGeometrySection(FEBioImport* pim) : FEBioFileSection(pim) {}

protected:
	bool ReadElement(XMLTag& tag, FEElement& el, int nid);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection1x : public FEBioGeometrySection
{
protected:
	struct FEDOMAIN
	{
		FE_Element_Spec		elem;	// element type
		int					mat;	// material ID
		int					nel;	// number of elements
	};

public:
	FEBioGeometrySection1x(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection(XMLTag& tag);
	void ParseNodeSetSection(XMLTag& tag);
	void ParseElementSection(XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection2 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection2(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection(XMLTag& tag);
	void ParseEdgeSection(XMLTag& tag);
	void ParseSurfaceSection   (XMLTag& tag);
	void ParseElementSection   (XMLTag& tag);
	void ParseNodeSetSection   (XMLTag& tag);
	void ParseElementSetSection(XMLTag& tag);
	void ParseElementDataSection(XMLTag& tag);
	void ParseElementData(FEElement& el, XMLTag& tag);
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection25 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection25(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseNodeSetSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseInstanceSection   (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);

	// New functions for parsing parts
	void ParsePart(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartElementSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSetSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartSurfaceSection(XMLTag& tag, FEBModel::Part* part);

protected:
	FEBModel			m_feb;
};

//-----------------------------------------------------------------------------
class FEBioGeometrySection3 : public FEBioGeometrySection
{
public:
	FEBioGeometrySection3(FEBioImport* pim) : FEBioGeometrySection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodeSection       (XMLTag& tag);
	void ParseDiscreteSetSection(XMLTag& tag);
	void ParseSurfacePairSection(XMLTag& tag);
	void ParseNodeSetPairSection(XMLTag& tag);
	void ParseNodeSetSetSection (XMLTag& tag);
	void ParsePartSection       (XMLTag& tag);
	void ParseInstanceSection   (XMLTag& tag);
	void ParseSurfaceSection    (XMLTag& tag);
	void ParseElementSection    (XMLTag& tag);
	void ParseNodeSetSection    (XMLTag& tag);
	void ParseEdgeSection       (XMLTag& tag);
	void ParseElementSetSection (XMLTag& tag);

	// New functions for parsing parts
	void ParsePart(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartElementSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartNodeSetSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartSurfaceSection(XMLTag& tag, FEBModel::Part* part);
	void ParsePartElementSetSection(XMLTag& tag, FEBModel::Part* part);

protected:
	FEBModel			m_feb;
};
