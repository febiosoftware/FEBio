#pragma once
#include <FECore/FEPointFunction.h>
#include <vector>
#include <string>
#include "FEDataSource.h"
using namespace std;

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
	virtual double Evaluate(vector<double>& f);

	// evaluate objective function
	double Evaluate();

	// print output to screen or not
	void SetVerbose(bool b) { m_verbose = b; }

	// return the FE model
	FEModel* GetFEModel() { return m_fem; }

public: // These functions need to be implemented by derived classes

	// return number of measurements (i.e. nr of terms in objective function)
	virtual int Measurements() = 0;

	// evaluate the function values (i.e. the f_i above)
	virtual void EvaluateFunctions(vector<double>& f) = 0;

	// get the measurement vector (i.e. the y_i above)
	virtual void GetMeasurements(vector<double>& y) = 0;

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
	void SetMeasurements(const vector<pair<double, double> >& data);

public:
	// return number of measurements
	int Measurements();

	// evaluate the function values
	void EvaluateFunctions(vector<double>& f);

	// get the measurement vector
	void GetMeasurements(vector<double>& y);

private:
	FEPointFunction		m_lc;		//!< data load curve for evaluating measurements
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
	void EvaluateFunctions(vector<double>& f) override;

	// get the measurement vector
	void GetMeasurements(vector<double>& y) override;

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

	void SetVariable(int n);

public:
	// return number of measurements (i.e. nr of terms in objective function)
	int Measurements() override;

	// evaluate the function values (i.e. the f_i above)
	void EvaluateFunctions(vector<double>& f) override;

	// get the measurement vector (i.e. the y_i above)
	void GetMeasurements(vector<double>& y) override;

private:
	std::vector<Entry>	m_Data;
	int					m_var;
};
