#pragma once
#include <FECore/FEPointFunction.h>

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

	// Evaluate source at time t
	virtual double Evaluate(double t) = 0;

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

	// Initialize data
	bool Init() override;

	// Reset data
	void Reset() override;

	// Evaluate the model parameter at time t
	double Evaluate(double t) override;

private:
	static bool update(FEModel* pmdl, unsigned int nwhen, void* pd);
	void update();

private:
	string	m_name;				//!< name of parameter that generates the function data
	double*	m_pd;				//!< pointer to variable data
	FEPointFunction		m_rf;	//!< reaction force data
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

	// evaluate data source at time t
	double Evaluate(double t) override;

private:
	FEDataSource*	m_src;
};
