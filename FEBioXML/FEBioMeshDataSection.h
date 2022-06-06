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
class FEDomainMap;
class FESurfaceMap;
class FENodeDataMap;

//-----------------------------------------------------------------------------
// MeshData Section (introduced in febio_spec 2.5)
class FEBioMeshDataSection : public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[FEElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection(FEBioImport* pim) : FEBioFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	void ParseShellThickness(XMLTag& tag, FEElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, FEElementSet& set);
	void ParseMaterialAxes  (XMLTag& tag, FEElementSet& set);
	void ParseMaterialData  (XMLTag& tag, FEElementSet& set, const string& name);
	void ParseMaterialFiberProperty(XMLTag& tag, FEElementSet& set);

private:
	void ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseDataArray(XMLTag& tag, FEDataArray& map, const char* sztag);
	FEElement* ParseElement(XMLTag& tag, FEElementSet& set);

private:
	vector<FEElement*> m_pelem;
};

//-----------------------------------------------------------------------------
// MeshData Section for febio_spec 3.0
class FEBioMeshDataSection3 : public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[FEElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection3(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalData(XMLTag& tag);
	void ParseSurfaceData(XMLTag& tag);
	void ParseElementData(XMLTag& tag);

protected:
//	void ParseModelParameter(XMLTag& tag, FEParamValue param);
//	void ParseMaterialPointData(XMLTag& tag, FEParamValue param);

protected:
	void ParseShellThickness(XMLTag& tag, FEElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, FEElementSet& set);
	void ParseMaterialAxes(XMLTag& tag, FEElementSet& set);
	void ParseMaterialAxesProperty(XMLTag& tag, FEElementSet& set);

private:
	void ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseSurfaceData(XMLTag& tag, FESurfaceMap& map);
	void ParseNodeData   (XMLTag& tag, FENodeDataMap& map);
};

//-----------------------------------------------------------------------------
// MeshData Section for febio_spec 4.0
class FEBioMeshDataSection4: public FEBioFileSection
{
	struct ELEMENT_DATA
	{
		int		nval;	// number of values read
		double	val[FEElement::MAX_NODES];	// scalar value
	};

public:
	FEBioMeshDataSection4(FEBioImport* pim) : FEBioFileSection(pim) {}
	void Parse(XMLTag& tag);

protected:
	void ParseNodalData(XMLTag& tag);
	void ParseSurfaceData(XMLTag& tag);
	void ParseElementData(XMLTag& tag);

private:
	void ParseNodeData(XMLTag& tag, FENodeDataMap& map);
	void ParseSurfaceData(XMLTag& tag, FESurfaceMap& map);
	void ParseElementData(XMLTag& tag, FEElementSet& set, vector<ELEMENT_DATA>& values, int nvalues);
	void ParseElementData(XMLTag& tag, FEDomainMap& map);
	void ParseShellThickness(XMLTag& tag, FEElementSet& set);
	void ParseMaterialAxes(XMLTag& tag, FEElementSet& set);
	void ParseMaterialFibers(XMLTag& tag, FEElementSet& set);
};
