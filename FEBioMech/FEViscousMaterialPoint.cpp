/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, Columbia University, and others.

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
#include "FEViscousMaterialPoint.h"
#include "FEElasticMaterial.h"
#include <FECore/DumpStream.h>

FEMaterialPoint* FEViscousMaterialPoint::Copy()
{
    FEViscousMaterialPoint* pt = new FEViscousMaterialPoint(*this);
    if (m_pNext) pt->m_pNext = m_pNext->Copy();
    return pt;
}

void FEViscousMaterialPoint::Init()
{
	// initialize data to identity
    m_Fp = mat3dd(1);

	m_dt = 0.0;
    
    // don't forget to intialize the nested data
    FEMaterialPoint::Init();
}

void FEViscousMaterialPoint::Update(const FETimeInfo& timeInfo)
{
	FEElasticMaterialPoint& pt = *this->ExtractData<FEElasticMaterialPoint>();
    m_Fp = pt.m_F;

	m_dt = timeInfo.timeIncrement;
    
    // don't forget to call base class
    FEMaterialPoint::Update(timeInfo);
}

void FEViscousMaterialPoint::Serialize(DumpStream& ar)
{
	FEMaterialPoint::Serialize(ar);
	ar & m_Fp;
}

mat3d FEViscousMaterialPoint::VelocityGradient()
{
    FEElasticMaterialPoint& pt = *this->ExtractData<FEElasticMaterialPoint>();
    mat3d Fi = pt.m_F.inverse();
    mat3d L = (mat3dd(1) - m_Fp*Fi)/m_dt;
    return L;
}

mat3ds FEViscousMaterialPoint::RateOfDeformation()
{
    mat3d L = VelocityGradient();
    return L.sym();
}
