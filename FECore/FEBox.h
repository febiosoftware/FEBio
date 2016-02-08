#pragma once
#include "vec3d.h"

//-----------------------------------------------------------------------------
class FEMesh;
class FEDomain;

//-----------------------------------------------------------------------------
// Helper class for finding bounding boxes.
class FEBox
{
public:
	FEBox();
	FEBox(const vec3d& r0, const vec3d& r1);

	FEBox(const FEMesh& mesh);
	FEBox(const FEDomain& dom);

	vec3d center() { return (m_r0 + m_r1)*0.5; }
	
	double width () { return m_r1.x - m_r0.x; };	//!< x-size
	double height() { return m_r1.y - m_r0.y; };	//!< y-size
	double depth () { return m_r1.z - m_r0.z; };	//!< z-size

	// the maximum size
	double maxsize();

	// union of box and node
	void add(const vec3d& r);

private:
	vec3d	m_r0, m_r1;
};
