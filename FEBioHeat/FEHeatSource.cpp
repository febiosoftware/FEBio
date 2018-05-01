#include "FEHeatSource.h"

//-----------------------------------------------------------------------------
FEHeatSource::FEHeatSource(FEModel* pfem) : FEBodyLoad(pfem)
{
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEConstHeatSource, FEHeatSource);
	ADD_PARAMETER(m_Q, FE_PARAM_DOUBLE, "Q");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEConstHeatSource::FEConstHeatSource(FEModel* pfem) : FEHeatSource(pfem)
{
	m_Q = 0;
}
