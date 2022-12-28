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
#include <FECore/PointCurve.h>
#include <functional>
#include <FECore/NodeDataRecord.h>

//-------------------------------------------------------------------------------------------------
// The FEDataSource class is used by the FEObjectiveFunction to query model data and evaluate it
// at the requested time point. This is an abstract base class and derived classes must implement
// the Evaluate function. 
class FEDataSource
{
public:
	FEDataSource(FEModel* fem);
	virtual ~FEDataSource();

	// Initialize data source
	virtual bool Init();

	// Reset data source 
	virtual void Reset();

	// Evaluate source at x
	virtual double Evaluate(double x) = 0;

protected:
	FEModel&			m_fem;	//!< reference to model
};

//-------------------------------------------------------------------------------------------------
// The FEDataParameter class is a data source the extracts data from a model parameter. The parameter
// must be set with SetParameterName before calling Init. 
class FEDataParameter : public FEDataSource
{
public:
	// constructor
	FEDataParameter(FEModel* fem);
	
	// Set the model parameter name
	void SetParameterName(const std::string& name);

	// set the ordinate name
	void SetOrdinateName(const std::string& name);

	// Initialize data
	bool Init() override;

	// Reset data
	void Reset() override;

	// Evaluate the model parameter at x
	double Evaluate(double x) override;

	// evaluate the current value
	double value() { return m_fy(); }

private:
	static bool update(FEModel* pmdl, unsigned int nwhen, void* pd);
	void update();

private:
	string	m_param;			//!< name of parameter that generates the function data
	string	m_ord;				//!< name of ordinate parameter
	std::function<double()>	m_fx;				//!< pointer to ordinate value
	std::function<double()>	m_fy;				//!< pointer to variable data
	PointCurve		m_rf;	//!< reaction force data
};

//-------------------------------------------------------------------------------------------------
// This data source class evaluates the data using another data source and then applying a filter. 
// In this case, the filter only returns positive values or zero otherwise. The data source must be
// set before calling Init
class FEDataFilterPositive : public FEDataSource
{
public:
	FEDataFilterPositive(FEModel* fem);
	~FEDataFilterPositive();

	// Set the data source
	void SetDataSource(FEDataSource* src);

	// Initialize data
	bool Init() override;

	// reset data
	void Reset() override;

	// evaluate data source at x
	double Evaluate(double x) override;

private:
	FEDataSource*	m_src;
};

//-------------------------------------------------------------------------------------------------
class FEDataFilterSum : public FEDataSource
{
public:
	FEDataFilterSum(FEModel* fem);
	~FEDataFilterSum();

	void SetData(FELogNodeData* data, FENodeSet* nodeSet);

	// Initialize data
	bool Init() override;

	// reset data
	void Reset() override;

	// evaluate data source at x
	double Evaluate(double x) override;


private:
	static bool update(FEModel* pmdl, unsigned int nwhen, void* pd);
	void update();

private:
	FELogNodeData*	m_data;
	FENodeSet*		m_nodeSet;
	PointCurve		m_rf;
};
