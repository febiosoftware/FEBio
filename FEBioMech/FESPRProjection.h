#pragma once
#include <vector>

class FESolidDomain;

//-------------------------------------------------------------------------------------------------
//! This class implements the super-convergent-patch recovery method which projects integration point
//! data to the finite element nodes.
class FESPRProjection
{
public:
	FESPRProjection();

	void Project(FESolidDomain& dom, const std::vector< std::vector<double> >& d, std::vector<double>& o);

	void SetInterpolationOrder(int p);

protected:
	int		m_p;	//!< interpolation order (set to -1 for default rules)
};
