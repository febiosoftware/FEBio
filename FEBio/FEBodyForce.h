#pragma once

//-----------------------------------------------------------------------------
//! This class is the base class for body forces
//
class FEBodyForce
{
public:
	FEBodyForce()
	{
		s[0] = s[1] = s[2] = 0.0;
		lc[0] = -1; lc[1] = -1; lc[2] = -1;
	}

public:
	double	s[3];		// scale factor
	int		lc[3];		// load curve number
};
