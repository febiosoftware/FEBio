/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2021 University of Utah, The Trustees of Columbia University in
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#pragma once

#include "mat3d.h"
#include "quatd.h"
#include "FETimeInfo.h"
#include <vector>

class FEElement;
class FEMaterialPoint;

//-----------------------------------------------------------------------------
//! Material point class

//! This class implements the concept of a material point. This point carries
//! with it not only information about its location, both in the reference and  
//! current configuration but also about the local deformation. In addition
//! it contains the state information that is associated with the current
//! point.
//! 

class FECORE_API FEMaterialPointData
{
public:
	FEMaterialPointData(FEMaterialPointData* ppt = 0);
	virtual ~FEMaterialPointData();

public:
	//! The init function is used to intialize data
	virtual void Init();

	//! The Update function is used to update material point data
	//! Note that this gets called at the start of the time step during PreSolveUpdate
	virtual void Update(const FETimeInfo& timeInfo);

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual FEMaterialPointData* Copy() { return 0; }

	// serialization
	virtual void Serialize(DumpStream& ar);

public:
	//! Get the next material point data
	FEMaterialPointData* Next() { return m_pNext; }

	//! Get the previous (parent) material point data
	FEMaterialPointData* Prev() { return m_pPrev; }
    
	// assign the previous pointer
	void SetPrev(FEMaterialPointData* pt);

	//! assign the next pointer
	//! this also sets the prev pointer of the passed pointer
	//! in other words, it makes this the parent of the passed pointer
	void SetNext(FEMaterialPointData* pt);

	//! append a material point
	void Append(FEMaterialPointData* pt);

protected:
	virtual int Components() const { return 0; }
	virtual FEMaterialPoint* GetPointData(int i) { return nullptr; }

public:
	//! Extract data (\todo Is it safe for a plugin to use this function?)
	template <class T> T* ExtractData();
	template <class T> const T* ExtractData() const;

protected:
	FEMaterialPointData*	m_pNext;    //!< next data in the list
	FEMaterialPointData*	m_pPrev;    //!< previous data in the list

	friend class FEMaterialPoint;
};

//-----------------------------------------------------------------------------
class FECORE_API FEMaterialPoint
{
public:
	FEMaterialPoint(FEMaterialPointData* data = nullptr);
	virtual ~FEMaterialPoint();

	//! The init function is used to intialize data
	virtual void Init();

	virtual FEMaterialPoint* Copy();

	//! The Update function is used to update material point data
	//! Note that this gets called at the start of the time step during PreSolveUpdate
	virtual void Update(const FETimeInfo& timeInfo);

	virtual void Serialize(DumpStream& ar);

	void Append(FEMaterialPointData* pt);

public:
	int Components() const { return (m_data ? m_data->Components() : 0); }

	FEMaterialPoint* GetPointData(int i = 0)
	{
		if (m_data == nullptr) return this;
		FEMaterialPoint* mp = m_data->GetPointData(i);
		if (mp == nullptr) mp = this;
		return mp;
	}

public:
	//! Extract data (\todo Is it safe for a plugin to use this function?)
	template <class T> T* ExtractData();
	template <class T> const T* ExtractData() const;

public:
	vec3d		m_r0;		//!< material point position
	vec3d		m_rt;		//!< current point position
	double		m_J0;		//!< reference Jacobian
	double		m_Jt;		//!< current Jacobian
	quatd		m_Q;		//!< local coordinates
	FEElement* m_elem;		//!< Element where this material point is
	int			m_index;	//!< local integration point index 

	// pointer to element's shape function values
	double* m_shape;

protected:
	FEMaterialPointData* m_data;
};

//-----------------------------------------------------------------------------
template <class T> inline T* FEMaterialPointData::ExtractData()
{
	// first see if this is the correct type
	T* p = dynamic_cast<T*>(this);
	if (p) return p;

	// check all the child classes 
	FEMaterialPointData* pt = this;
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

	if (Components() > 0)
	{
		for (int i = 0; i < Components(); ++i)
		{
			FEMaterialPoint* mpi = GetPointData(i);
			p = mpi->ExtractData<T>();
			if (p) return p;
		}
	}

	// Everything has failed. Material point data can not be found
	return 0;
}

//-----------------------------------------------------------------------------
template <class T> inline const T* FEMaterialPointData::ExtractData() const
{
	// first see if this is the correct type
	const T* p = dynamic_cast<const T*>(this);
	if (p) return p;

	// check all the child classes 
	const FEMaterialPointData* pt = this;
	while (pt->m_pNext)
	{
		pt = pt->m_pNext;
		p = dynamic_cast<const T*>(pt);
		if (p) return p;
	}

	// search up
	pt = this;
	while (pt->m_pPrev)
	{
		pt = pt->m_pPrev;
		p = dynamic_cast<const T*>(pt);
		if (p) return p;
	}

	// Everything has failed. Material point data can not be found
	return 0;
}

//-----------------------------------------------------------------------------
template <class T> inline T* FEMaterialPoint::ExtractData()
{
	return (m_data ? m_data->ExtractData<T>() : nullptr);
}

//-----------------------------------------------------------------------------
template <class T> inline const T* FEMaterialPoint::ExtractData() const
{
	return (m_data ? m_data->ExtractData<T>() : nullptr);
}

//-----------------------------------------------------------------------------
// Material point base class for materials that define vector properties
class FECORE_API FEMaterialPointArray : public FEMaterialPointData
{
public:
	FEMaterialPointArray(FEMaterialPointData* ppt = nullptr);

	//! initialization
	void Init() override;

	//! serialization
	void Serialize(DumpStream& ar) override;

	//! material point update
	void Update(const FETimeInfo& timeInfo) override;

public:
	//! Add a child material point
	void AddMaterialPoint(FEMaterialPoint* pt);

	//! get the number of material point components
	int Components() const override { return (int)m_mp.size(); }

	//! retrieve point data
	FEMaterialPoint* GetPointData(int i) override { return m_mp[i]; }

protected:
	std::vector<FEMaterialPoint*>	m_mp;	//!< material point data for indidivual properties
};
