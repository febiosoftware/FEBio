#pragma once

#include "FECore/FEAnalysis.h"

class FEModel;

//-----------------------------------------------------------------------------
class FEBioProgress : public Progress
{
public:
	FEBioProgress(FEModel& fem) : m_fem(fem) {}
	void SetProgress(double f);

protected:
	FEModel&	m_fem;
};
