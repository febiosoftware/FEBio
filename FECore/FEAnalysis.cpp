#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"

FEDomain* FEAnalysis::Domain(int i)
{
	return &(m_fem.m_mesh.Domain(m_Dom[i])); 
}
