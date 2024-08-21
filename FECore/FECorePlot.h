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
#include "FEPlotData.h"
#include "fecore_api.h"

//-----------------------------------------------------------------------------
class FEMaterial;
class FEDomainList;
class FEFacetSet;
class FEDataMap;

//-----------------------------------------------------------------------------
// class for exporting element specific material parameters to plot file
class FECORE_API FEPlotParameter : public FEPlotData
{
	FECORE_BASE_CLASS(FEPlotParameter)

public:
	FEPlotParameter(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a) override;
	bool Save(FESurface& dom, FEDataStream& a) override;
	bool Save(FEMesh& mesh, FEDataStream& a) override;

	bool SetFilter(const char* sz) override;

	void Serialize(DumpStream& ar) override;

protected:
	FEParamValue	m_param;	//!< parameter
	int				m_index;	//!< index for array parameters

private:
	FEMaterial*			m_mat;
	FEDomainList*		m_dom;
	FEFacetSet*			m_surf;

	std::string		m_filter;
};

//-----------------------------------------------------------------------------
class FEPIDController;

class FECORE_API FEPlotPIDController : public FEPlotGlobalData
{
public:
	FEPlotPIDController(FEModel* pfem);
	bool SetFilter(const char* sz) override;
	bool Save(FEDataStream& a) override;

private:
	FEPIDController* m_pid;
};

//-----------------------------------------------------------------------------
class FECORE_API FEPlotMeshData : public FEPlotData
{
	FECORE_BASE_CLASS(FEPlotMeshData)

public:
	FEPlotMeshData(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a) override;
	bool Save(FESurface& dom, FEDataStream& a) override;
	bool Save(FEMesh& mesh, FEDataStream& a) override;

	bool SetFilter(const char* sz) override;

	void Serialize(DumpStream& ar) override;

protected:
	FEDataMap*	m_map;		//!< parameter

	const FEDomainList*	m_dom;
	std::string		m_filter;
};

//! Plot variable for field variables
class FECORE_API FEPlotFieldVariable : public FEPlotNodeData
{
public:
	FEPlotFieldVariable(FEModel* pfem);
	bool Save(FEMesh& mesh, FEDataStream& a) override;

	bool SetFilter(const char* sz) override;

	void Serialize(DumpStream& ar) override;

protected:
	std::vector<int>	m_dofs;
};
