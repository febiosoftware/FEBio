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
#include "FECoreBase.h"
#include "FEMaterialPoint.h"
#include "fecore_type.h"

//--------------------------------------------------------------------
// This class defines a material point property. Material point properties
// are used to move data in and out of material points.
class FEMaterialPointProperty
{
public:
	FEMaterialPointProperty(FEDataType dataType) : m_type(dataType) {}
	virtual ~FEMaterialPointProperty() {}

	FEDataType dataType() const { return m_type; }

	virtual void set(FEMaterialPoint& mp, const double& v) {}
	virtual void set(FEMaterialPoint& mp, const vec3d&  v) {}
	virtual void set(FEMaterialPoint& mp, const mat3d&  v) {}

	virtual void get(FEMaterialPoint& mp, double& v) {}
	virtual void get(FEMaterialPoint& mp, vec3d&  v) {}
	virtual void get(FEMaterialPoint& mp, mat3d&  v) {}

private:
	FEDataType	m_type;
};

//--------------------------------------------------------------------
// Template class for defining material point property classes. 
// Derived classes need to override the data member. 
template <class T, class A> class FEMaterialPointProperty_T : public FEMaterialPointProperty
{
public:
	FEMaterialPointProperty_T() : FEMaterialPointProperty(fecoreType<A>::type()) {}
	void set(FEMaterialPoint& mp, const A& Q)
	{
		T& pt = *mp.ExtractData<T>();
		data(pt) = Q;
	}

	void get(FEMaterialPoint& mp, A& Q)
	{
		T& pt = *mp.ExtractData<T>();
		Q = data(pt);
	}

	virtual A& data(T& pt) = 0;
};
