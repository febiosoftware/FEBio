#pragma once

#include "mat3d.h"
#include "FEParameterList.h"
#include "FETypes.h"
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
	virtual void Init();

	//! The Update function is used to update material point data
	//! Note that this gets called at the start of the time step during PreSolveUpdate
	virtual void Update(const FETimeInfo& timeInfo);

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual FEMaterialPoint* Copy() = 0;

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

	//! assign the next pointer
	//! this also sets the prev pointer of the passed pointer
	//! in other words, it makes this the parent of the passed pointer
	void SetNext(FEMaterialPoint* pt);

	// serialization
	void Serialize(DumpStream& ar);

	// find a parameter with a given name
	FEParam* FindParameter(const char* szname);

protected:
	FEMaterialPoint*	m_pNext;	//<! next data in the list
	FEMaterialPoint*	m_pPrev;	//<! previous data in the list
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
