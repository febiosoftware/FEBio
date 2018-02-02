#pragma once
#include <FECore/FEDataLoadCurve.h>

class FEDataSource
{
public:
	FEDataSource(FEModel* fem);
	virtual ~FEDataSource();

	virtual bool Init();

	virtual void Reset();

	virtual double Evaluate(double t) = 0;

protected:
	FEModel&			m_fem;
};

class FEDataParameter : public FEDataSource
{
public:
	FEDataParameter(FEModel* fem);

	void SetParameterName(const std::string& name);

	bool Init() override;

	void Reset() override;

	double Evaluate(double t) override;

private:
	static bool update(FEModel* pmdl, unsigned int nwhen, void* pd);
	void update();

private:
	string	m_name;				//!< name of parameter that generates the function data
	double*	m_pd;				//!< pointer to variable data
	FEDataLoadCurve		m_rf;	//!< reaction force data
};

class FEDataFilterPositive : public FEDataSource
{
public:
	FEDataFilterPositive(FEModel* fem);
	~FEDataFilterPositive();

	void SetDataSource(FEDataSource* src);

	bool Init() override;

	void Reset() override;

	double Evaluate(double t) override;

private:
	FEDataSource*	m_src;
};
