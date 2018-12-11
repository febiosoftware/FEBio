#include "stdafx.h"
#include "FEDataSource.h"
#include <FECore/FEModel.h>
#include <FECore/fecore_error.h>

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
	double x = *(m_px);
	double y = *(m_py);

	// add the data pair to the loadcurve
	m_rf.Add(x, y);
}

FEDataParameter::FEDataParameter(FEModel* fem) : FEDataSource(fem), m_rf(fem)
{
	m_ord = "fem.time";
}

void FEDataParameter::SetParameterName(const std::string& name)
{
	m_param = name;
}

// set the ordinate name
void FEDataParameter::SetOrdinateName(const std::string& name)
{
	m_ord = name;
}

bool FEDataParameter::Init()
{
	// find all the parameters
	FEParamValue val = m_fem.GetParameterValue(ParamString(m_param.c_str()));
	if (val.isValid() == false) return fecore_error("Invalid parameter name %s", m_param.c_str());
	if (val.type() != FE_PARAM_DOUBLE) return fecore_error("Invalid type for parameter %s", m_param.c_str());
	m_py = (double*)val.data_ptr();
	if (m_py == 0) return fecore_error("Invalid data pointer for parameter %s", m_param.c_str());

	// find the ordinate
	val = m_fem.GetParameterValue(ParamString(m_ord.c_str()));
	if (val.isValid() == false) return fecore_error("Invalid ordinate name %s", m_ord.c_str());
	if (val.type() != FE_PARAM_DOUBLE) return fecore_error("Invalid type for ordinate %s", m_ord.c_str());
	m_px = (double*)val.data_ptr();
	if (m_px == 0) return fecore_error("Invalid data pointer for ordinate %s", m_ord.c_str());

	// register callback
	m_fem.AddCallback(update, CB_INIT | CB_MAJOR_ITERS, (void*) this);

	return FEDataSource::Init();
}

void FEDataParameter::Reset()
{
	// reset the reaction force load curve
	m_rf.Clear();
	*m_px = 0.0;
	FEDataSource::Reset();
}

double FEDataParameter::Evaluate(double x)
{
	return m_rf.value(x);
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
