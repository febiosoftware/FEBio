#pragma once

#include "mat3d.h"
#include "DumpFile.h"
#include "DumpStream.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Material point class

//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and  
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.

class FEMaterialPoint
{
public:
	FEMaterialPoint(FEMaterialPoint* ppt = 0) : m_pt(ppt){}
	virtual ~FEMaterialPoint() { if (m_pt) { delete m_pt; m_pt = 0; } }

	//! The init function is used to intialize data
	virtual void Init(bool bflag) = 0;

	virtual FEMaterialPoint* Copy() = 0;

	virtual void Serialize(DumpFile& ar) = 0;

	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	template <class T> T* ExtractData();

	virtual FEMaterialPoint* GetPointData(int i) { return this; }

protected:
	FEMaterialPoint*	m_pt;	//<! nested point data

public:
	static double time;	// time value
	static double dt; // time increment
};

//-----------------------------------------------------------------------------
template <class T> inline T* FEMaterialPoint::ExtractData()
{
	T* p = dynamic_cast<T*>(this);
	if (p) return p; 
	else
	{
		if (m_pt) return m_pt->ExtractData<T>();
		else return 0;
	}
}
