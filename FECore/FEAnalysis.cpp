#include "stdafx.h"
#include "FEAnalysis.h"
#include "FEModel.h"

FEDomain* FEAnalysis::Domain(int i)
{
	return &(m_fem.GetMesh().Domain(m_Dom[i])); 
}
