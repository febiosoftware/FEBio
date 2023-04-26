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



#include "stdafx.h"
#include "FEMat3dsValuator.h"
#include "FEModel.h"
#include "FEMesh.h"
#include "FEDataMap.h"
#include "DumpStream.h"

//=============================================================================
// FEConstValueMat3ds
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEConstValueMat3ds, FEMat3dsValuator)
	ADD_PARAMETER(m_val, "const");
END_FECORE_CLASS();

FEConstValueMat3ds::FEConstValueMat3ds(FEModel* fem) : FEMat3dsValuator(fem)
{
	m_val.zero();
}

FEMat3dsValuator* FEConstValueMat3ds::copy()
{
	FEConstValueMat3ds* map = fecore_alloc(FEConstValueMat3ds, GetFEModel());
	map->m_val = m_val;
	return map;
}

//=============================================================================
// FEMappedValueMat3ds
//-----------------------------------------------------------------------------

BEGIN_FECORE_CLASS(FEMappedValueMat3ds, FEMat3dsValuator)
	ADD_PARAMETER(m_mapName, "map");
END_FECORE_CLASS();

FEMappedValueMat3ds::FEMappedValueMat3ds(FEModel* fem) : FEMat3dsValuator(fem)
{
	m_val = nullptr;
}

void FEMappedValueMat3ds::setDataMap(FEDataMap* val)
{
	m_val = val;
}

mat3ds FEMappedValueMat3ds::operator()(const FEMaterialPoint& pt)
{
	return m_val->valueMat3ds(pt);
}

FEMat3dsValuator* FEMappedValueMat3ds::copy()
{
	FEMappedValueMat3ds* map = fecore_alloc(FEMappedValueMat3ds, GetFEModel());
	map->m_val = m_val;
	return map;
}

void FEMappedValueMat3ds::Serialize(DumpStream& ar)
{
	FEMat3dsValuator::Serialize(ar);
	ar & m_val;
}
