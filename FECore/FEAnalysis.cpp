#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"

FEAnalysis::FEAnalysis(FEModel& fem, int ntype) : m_fem(fem), m_ntype(ntype) 
{
	m_tend = 0.0;
}

FEDomain* FEAnalysis::Domain(int i)
{
	return &(m_fem.GetMesh().Domain(m_Dom[i])); 
}
