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
	FEMaterialPoint(FEMaterialPoint* ppt = 0);
	virtual ~FEMaterialPoint();

public:
	//! The init function is used to intialize data
	virtual void Init(bool bflag);

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual FEMaterialPoint* Copy() = 0;

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual void ShallowCopy(DumpStream& dmp, bool bsave) = 0;

	//! get the number of material point components
	virtual int Components() { return 1; }

	//! Get the material point data
	virtual FEMaterialPoint* GetPointData(int i) { return this; }

	//! Get the next material point data
	FEMaterialPoint* Next() { return m_pNext; }

	//! Get the previous (parent) material point data
	FEMaterialPoint* Prev() { return m_pPrev; }
    
	//! Extract data (\todo Is it safe for a plugin to use this function?)
	template <class T> T* ExtractData();

	// assign the previous pointer
	void SetPrev(FEMaterialPoint* pt);

protected:
	FEMaterialPoint*	m_pNext;	//<! next data in the list
	FEMaterialPoint*	m_pPrev;	//<! previous data in the list

public:
	static double time;	// time value
	static double dt; // time increment
};

//-----------------------------------------------------------------------------
template <class T> inline T* FEMaterialPoint::ExtractData()
{
	// first see if this is the correct type
	T* p = dynamic_cast<T*>(this);
	if (p) return p;

	// check all the child classes 
	FEMaterialPoint* pt = this;
	while (pt->m_pNext)
	{
		pt = pt->m_pNext;
		p = dynamic_cast<T*>(pt);
		if (p) return p;
	}

	// search up
	pt = this;
	while (pt->m_pPrev)
	{
		pt = pt->m_pPrev;
		p = dynamic_cast<T*>(pt);
		if (p) return p;
	}

	// Everything has failed. Material point data can not be found
	return 0;
}
