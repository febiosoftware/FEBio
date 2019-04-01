#pragma once
#include "vec3d.h"
#include <vector>
#include "fecore_api.h"

class FESurface;

//-----------------------------------------------------------------------------
//! This class is a helper class to locate the nearest neighbour on a surface

class FECORE_API FENNQuery
{
public:
	struct NODE
	{
		int		i;	// index of node
		vec3d	r;	// position of node
		double	d1;	// distance to pivot 1
		double	d2;	// distance to pivot 2
	};

public:
	FENNQuery(FESurface* ps = 0);
	virtual ~FENNQuery();

	//! initialize search structures
	void Init();
	void InitReference();

	//! attach to a surface
	void Attach(FESurface* ps) { m_ps = ps; }

	//! find the neirest neighbour of r
	int Find(vec3d x);	
	int FindReference(vec3d x);	

protected:
	int FindRadius(double r);

protected:
	FESurface*	m_ps;	//!< the surface to search
	std::vector<NODE>	m_bk;	// BK tree

	vec3d	m_q1;	// pivot 1
	vec3d	m_q2;	// pivot 2

	int		m_imin;	// last found index
};
