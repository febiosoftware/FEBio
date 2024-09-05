#pragma once
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
#include <FECore/FEPlotData.h>
#include <FECore/FEElement.h>
#include <FECore/units.h>

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! growth tensor average value
class FEPlotGrowthRatio : public FEPlotDomainData
{
public:
	FEPlotGrowthRatio(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! growth tensor average value
class FEPlotGrowthTensor : public FEPlotDomainData
{
public:
	FEPlotGrowthTensor(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Stress traces
class FEPlotTraceStresses : public FEPlotDomainData
{
public:
	FEPlotTraceStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fe
class FEPlotGrowthElasticDeformationGradient : public FEPlotDomainData
{
public:
	FEPlotGrowthElasticDeformationGradient(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Fg=FFe^-1
class FEPlotGrowthDeformationGradient : public FEPlotDomainData
{
public:
	FEPlotGrowthDeformationGradient(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! J=det(F)
class FEPlotJacobian : public FEPlotDomainData
{
public:
	FEPlotJacobian(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Je=det(Fe)
class FEPlotGrowthElasticJacobian : public FEPlotDomainData
{
public:
	FEPlotGrowthElasticJacobian(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Jg=det(Fg)
class FEPlotGrowthJacobian : public FEPlotDomainData
{
public:
	FEPlotGrowthJacobian(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Plot k(theta)
class FEPlotGrowthK : public FEPlotDomainData
{
public:
	FEPlotGrowthK(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Growth environmental function value (phi)
class FEPlotGrowthPhi : public FEPlotDomainData
{
public:
	FEPlotGrowthPhi(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Actual solute concentration (for elastic-reaction-diffusion problems)
class FEPlotActualSoluteConcentrationERD : public FEPlotDomainData
{
public:
	FEPlotActualSoluteConcentrationERD(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);
protected:
	vector<int>	m_sol;
};

//-----------------------------------------------------------------------------
//! Nodal effective solute concentrations (for elastic-reaction-diffusion problems)
class FEPlotEffectiveSoluteConcentrationERD : public FEPlotDomainData
{
public:
	FEPlotEffectiveSoluteConcentrationERD(FEModel* pfem);
	bool Save(FEDomain& m, FEDataStream& a);
protected:
	vector<int>	m_sol;
};

//-----------------------------------------------------------------------------
//! Solute flux (for elastic-reaction-diffusion problems)
class FEPlotSoluteFluxERD : public FEPlotDomainData
{
public:
	FEPlotSoluteFluxERD(FEModel* pfem);
	bool Save(FEDomain& dom, FEDataStream& a);

protected:
	vector<int>	m_sol;
};