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
#include <vector>
#include <string>
#include "FEDataSource.h"
#include <FECore/ElementDataRecord.h>
#include <FECore/NodeDataRecord.h>

class FEModel;
class FEElement;

//=============================================================================
//! This class evaluates the objective function, which is defined as the sum
//! of squares of non-linear functions. 
//            ----
//             \                   2
//  f_obj =    /    [ f_i(a) - y_i]
//            ----
//              i
//
//! It is an abstract base class and derived classes
//! need to implement the nonlinear functions.
class FEObjectiveFunction
{
public:
	// constructor
	FEObjectiveFunction(FEModel* fem);
	virtual ~FEObjectiveFunction();

	// one-time initialization
	virtual bool Init();

	// This is called before each optimization iteration
	// and should be used by derived classes to reset any data
	virtual void Reset();

	// evaluate objective function
	// also returns the function values in f
	virtual double Evaluate(std::vector<double>& f);

	// evaluate objective function
	double Evaluate();

	// print output to screen or not
	void SetVerbose(bool b) { m_verbose = b; }

	// return the FE model
	FEModel* GetFEModel() { return m_fem; }

public:
	double RegressionCoefficient(const std::vector<double>& y0, const std::vector<double>& y);

public: // These functions need to be implemented by derived classes

	// return number of measurements (i.e. nr of terms in objective function)
	virtual int Measurements() = 0;

	// evaluate the function values (i.e. the f_i above)
	virtual void EvaluateFunctions(std::vector<double>& f) = 0;

	// get the measurement vector (i.e. the y_i above)
	virtual void GetMeasurements(std::vector<double>& y) = 0;

private:
	FEModel*	m_fem;
	bool	m_verbose;		//!< print data flag
};

//=============================================================================
// Objective function for fitting value pairs data
// In this case the nonlinear functions are defined as:
// f_i(a) = F(t_i; a)
// where t_i is the time for time step i, and F is the function that will be fitted to the measurement vector.
class FEDataFitObjective : public FEObjectiveFunction
{
public:
	FEDataFitObjective(FEModel* fem);
	~FEDataFitObjective();

	// one-time initialization
	bool Init();

	// This is called before each optimization iteration
	// and should be used by derived classes to reset any data
	void Reset();

	// set the data source
	void SetDataSource(FEDataSource* src);

	// set the data measurements
	void SetMeasurements(const std::vector<pair<double, double> >& data);

public:
	// return number of measurements
	int Measurements();

	// evaluate the function values
	void EvaluateFunctions(std::vector<double>& f);

	// get the measurement vector
	void GetMeasurements(std::vector<double>& y);

private:
	PointCurve			m_lc;		//!< data load curve for evaluating measurements
	FEDataSource*		m_src;		//!< source for evaluating functions
};

//=============================================================================
// Objective function for minimization of model parameters
class FEMinimizeObjective : public FEObjectiveFunction
{
	class Function
	{
	public:
		string		name;
		double*		var;

		double	y0;		// target value (i.e. "measurment")

	public:
		Function() : var(0), y0(0.0) {}
		Function(const Function& f) { name = f.name; var = f.var; y0 = f.y0; }
		void operator = (const Function& f) { name = f.name; var = f.var; y0 = f.y0; }
	};

public:
	FEMinimizeObjective(FEModel* fem);

	// one-time initialization
	bool Init() override;

	bool AddFunction(const char* szname, double targetValue);

public:
	// return number of measurements
	int Measurements() override;

	// evaluate the function values
	void EvaluateFunctions(std::vector<double>& f) override;

	// get the measurement vector
	void GetMeasurements(std::vector<double>& y) override;

private:
	std::vector<Function>	m_Func;
};

//=============================================================================
// This objective function evaluates element values with a user-defined
// table of values. The table stores for each element the required value. 
class FEElementDataTable : public FEObjectiveFunction
{
	struct Entry {
		int			elemId;	// the ID of the element
		double		target;	// the target value
		FEElement*	pe;		// pointer to element
	};

public:
	FEElementDataTable(FEModel* fem);

	bool Init() override;

	void AddValue(int elemID, double v);

	void SetVariable(FELogElemData* var);

public:
	// return number of measurements (i.e. nr of terms in objective function)
	int Measurements() override;

	// evaluate the function values (i.e. the f_i above)
	void EvaluateFunctions(std::vector<double>& f) override;

	// get the measurement vector (i.e. the y_i above)
	void GetMeasurements(std::vector<double>& y) override;

private:
	std::vector<Entry>	m_Data;
	FELogElemData*		m_var;
};

//=============================================================================
// This objective function evaluates node values with a user-defined
// table of values. The table stores for each node the required value. 
class FENodeDataTable : public FEObjectiveFunction
{
	struct Entry {
		int			nodeId;	// the ID of the node
		double		target;	// the target value
		int			index;	// zero-based node index
		int			ivar;  // variable
	};

public:
	FENodeDataTable(FEModel* fem);

	bool Init() override;

	bool AddValue(int elemID, std::vector<double>& v);

	void AddVariable(FELogNodeData* var);

public:
	// return number of measurements (i.e. nr of terms in objective function)
	int Measurements() override;

	// evaluate the function values (i.e. the f_i above)
	void EvaluateFunctions(std::vector<double>& f) override;

	// get the measurement vector (i.e. the y_i above)
	void GetMeasurements(std::vector<double>& y) override;

private:
	std::vector<Entry>				m_Data;
	std::vector<FELogNodeData*>		m_var;
};
