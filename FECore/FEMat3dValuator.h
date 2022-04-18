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
#include "FEValuator.h"

class FEDataMap;

//---------------------------------------------------------------------------------------
// Base class for evaluating vec3d parameters
class FECORE_API FEMat3dValuator : public FEValuator
{
	FECORE_SUPER_CLASS(FEMAT3DVALUATOR_ID)
	FECORE_BASE_CLASS(FEMat3dValuator)

public:
	FEMat3dValuator(FEModel* fem) : FEValuator(fem) {};

	virtual mat3d operator()(const FEMaterialPoint& pt) = 0;

	virtual FEMat3dValuator* copy() = 0;

	virtual bool isConst() { return false; }

	virtual mat3d* constValue() { return nullptr; }
};

//-----------------------------------------------------------------------------
// A constant valuator
class FECORE_API FEConstValueMat3d : public FEMat3dValuator
{
public:
	FEConstValueMat3d(FEModel* fem);

	FEMat3dValuator* copy() override;

	mat3d operator()(const FEMaterialPoint& pt) override { return m_val; }

	// is this a const value
	bool isConst() override { return true; }

	// get the const value (returns 0 if param is not const)
	mat3d* constValue() override { return &m_val; }

	mat3d& value() { return m_val; }

private:
	mat3d	m_val;

	DECLARE_FECORE_CLASS();
};


//-----------------------------------------------------------------------------
//! This class generates a material axes based on the local element node numbering.
class FECORE_API FEMat3dLocalElementMap : public FEMat3dValuator
{
public:
	FEMat3dLocalElementMap(FEModel* pfem);

	bool Init() override;

	mat3d operator () (const FEMaterialPoint& mp) override;

	FEMat3dValuator* copy() override;

	void Serialize(DumpStream& ar) override;

public:
	int			m_n[3];	// local element nodes

protected:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
//! This class generates material axes based on a spherical map. 
class FECORE_API FEMat3dSphericalMap : public FEMat3dValuator
{
public:
	FEMat3dSphericalMap(FEModel* pfem);

	bool Init() override;

	void SetSphereCenter(const vec3d& c) { m_c = c; }

	void SetSphereVector(const vec3d& r) { m_r = r; }

	mat3d operator() (const FEMaterialPoint& mp) override;

	FEMat3dValuator* copy() override;

public:
	vec3d		m_c;	// center of map
	vec3d		m_r;	// vector for parallel transport

protected:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEMat3dCylindricalMap : public FEMat3dValuator
{
public:
	FEMat3dCylindricalMap(FEModel* pfem);

	bool Init() override;

	void SetCylinderCenter(vec3d c) { m_c = c; }

	void SetCylinderAxis(vec3d a) { m_a = a; m_a.unit(); }

	void SetCylinderRef(vec3d r) { m_r = r; m_r.unit(); }

	mat3d operator () (const FEMaterialPoint& mp) override;

	FEMat3dValuator* copy() override;

public:
	vec3d		m_c;	// center of map
	vec3d		m_a;	// axis
	vec3d		m_r;	// reference direction

protected:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEMat3dPolarMap : public FEMat3dValuator
{
public:
	FEMat3dPolarMap(FEModel* pfem);

	bool Init() override;

	void SetCenter(vec3d c) { m_c = c; }

	void SetAxis(vec3d a) { m_a = a; m_a.unit(); }

	void SetVector0(vec3d r) { m_d0 = r; m_d0.unit(); }
	void SetVector1(vec3d r) { m_d1 = r; m_d1.unit(); }

	void SetRadius0(double r) { m_R0 = r; }
	void SetRadius1(double r) { m_R1 = r; }

	mat3d operator () (const FEMaterialPoint& mp) override;

	FEMat3dValuator* copy() override;

public:
	vec3d		m_c;		// center of map
	vec3d		m_a;		// axis
	vec3d		m_d0, m_d1;	// reference direction
	double		m_R0, m_R1;

protected:
	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FECORE_API FEMat3dVectorMap : public FEMat3dValuator
{
public:
	FEMat3dVectorMap(FEModel* pfem);

	bool Init() override;

	void SetVectors(vec3d a, vec3d d);

	mat3d operator () (const FEMaterialPoint& mp) override;

	FEMat3dValuator* copy() override;

	void Serialize(DumpStream& ar) override;

public:
	vec3d	m_a, m_d;
	mat3d	m_Q;

	DECLARE_FECORE_CLASS();
};

//---------------------------------------------------------------------------------------
class FECORE_API FEMappedValueMat3d : public FEMat3dValuator
{
public:
	FEMappedValueMat3d(FEModel* fem);

	void setDataMap(FEDataMap* val);

	FEDataMap* dataMap();

	mat3d operator()(const FEMaterialPoint& pt) override;

	FEMat3dValuator* copy() override;

	void Serialize(DumpStream& ar) override;

private:
	std::string	m_mapName;

private:
	FEDataMap* m_val;

	DECLARE_FECORE_CLASS();
};
