#pragma once

#include "mat3d.h"
#include "DumpFile.h"
#include "DumpStream.h"
#include "FEParameterList.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
//! Material point class

//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and  
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.

class FEMaterialPoint : public FEParamContainer
{
public:
	FEMaterialPoint(FEMaterialPoint* ppt = 0) : m_pt(ppt){}
	virtual ~FEMaterialPoint() { if (m_pt) { delete m_pt; m_pt = 0; } }

	//! The init function is used to intialize data
	virtual void Init(bool bflag) = 0;

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual FEMaterialPoint* Copy() = 0;

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	//! Get the material point data
	virtual FEMaterialPoint* GetPointData(int i) { return this; }

	//! Get the nested material point data
	FEMaterialPoint* Next() { return m_pt; }

	//! Extract data (\todo Is it safe for a plugin to use this function?)
	template <class T> T* ExtractData();

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
