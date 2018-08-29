#include "stdafx.h"
#include "FEDataSource.h"
#include <FECore/FEModel.h>

//=================================================================================================
FEDataSource::FEDataSource(FEModel* fem) : m_fem(*fem)
{
	
}

FEDataSource::~FEDataSource()
{
	
}

bool FEDataSource::Init()
{
	return true;
}

void FEDataSource::Reset()
{
	
}

//=================================================================================================
bool FEDataParameter::update(FEModel* pmdl, unsigned int nwhen, void* pd)
{
	// get the optimizaton data
	FEDataParameter& src = *((FEDataParameter*)pd);
	src.update();

	return true;
}

void FEDataParameter::update()
{
	// get the current time value
	double time = m_fem.GetTime().currentTime;

	// evaluate the current reaction force value
	double value = *(m_pd);

	// add the data pair to the loadcurve
	m_rf.Add(time, value);
}

FEDataParameter::FEDataParameter(FEModel* fem) : FEDataSource(fem), m_rf(fem)
{
}

void FEDataParameter::SetParameterName(const std::string& name)
{
	m_name = name;
}

bool FEDataParameter::Init()
{
	// find all the parameters
	FEParamValue val = m_fem.GetParameterValue(ParamString(m_name.c_str()));
	if (val.isValid() == false) return false;
	if (val.type() != FE_PARAM_DOUBLE) return false;
	m_pd = (double*)val.data_ptr();
	if (m_pd == 0) return false;

	// register callback
	m_fem.AddCallback(update, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return FEDataSource::Init();
}

void FEDataParameter::Reset()
{
	// reset the reaction force load curve
	m_rf.Clear();

	FEDataSource::Reset();
}

double FEDataParameter::Evaluate(double t)
{
	return m_rf.Value(t);
}

//=================================================================================================
FEDataFilterPositive::FEDataFilterPositive(FEModel* fem) : FEDataSource(fem)
{
	m_src = 0;
}

FEDataFilterPositive::~FEDataFilterPositive()
{
	if (m_src) delete m_src;
}

void FEDataFilterPositive::SetDataSource(FEDataSource* src)
{
	if (m_src) delete m_src;
	m_src = src;
}

bool FEDataFilterPositive::Init()
{
	if (m_src == 0) return false;

	return m_src->Init();
}

void FEDataFilterPositive::Reset()
{
	if (m_src) m_src->Reset();
}

double FEDataFilterPositive::Evaluate(double t)
{
	double v = m_src->Evaluate(t);
	return (v >= 0.0 ? v : -v);
}
