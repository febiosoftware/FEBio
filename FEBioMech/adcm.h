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
#include <FECore/ad.h>
#include <FECore/ad2.h>

namespace ad {

	template <class T>
	double StrainEnergy(T* p, FEMaterialPoint& mp)
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		::mat3ds C = pt.RightCauchyGreen();
		auto W = std::bind(&T::StrainEnergy_AD, p, mp, std::placeholders::_1);
		return ad::Evaluate(W, C);
	}

	template <class T>
	::mat3ds PK2Stress(T* p, FEMaterialPoint& mp, ::mat3ds& C)
	{
		auto W = std::bind(&T::StrainEnergy_AD, p, mp, std::placeholders::_1);
		return ad::Derive(W, C) * 2.0;
	}

	template <class T>
	::mat3ds PK2Stress(T* p, FEMaterialPoint& mp)
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		::mat3ds C = pt.RightCauchyGreen();
		return PK2Stress(p, mp, C);
	}

	template <class T>
	::tens4ds Tangent(T* p, FEMaterialPoint& mp, ::mat3ds& C)
	{
		auto S = std::bind(&T::PK2Stress_AD, p, mp, std::placeholders::_1);
		return Derive(S, C) * 2.0;
	}

	template <class T>
	::tens4ds Tangent(T* p, FEMaterialPoint& mp)
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		::mat3ds C = pt.RightCauchyGreen();
		return Tangent(p, mp, C);
	}
}

namespace ad2 {
	template <class T>
	::mat3ds PK2Stress(T* p, FEMaterialPoint& mp)
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		::mat3ds C = pt.RightCauchyGreen();
		auto W = std::bind(&T::StrainEnergy_AD2, p, mp, std::placeholders::_1);
		return ad2::Derive(W, C) * 2.0;
	}

	template <class T>
	::tens4ds Tangent(T* p, FEMaterialPoint& mp)
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		::mat3ds C = pt.RightCauchyGreen();
		auto W = std::bind(&T::StrainEnergy_AD2, p, mp, std::placeholders::_1);
		return ad2::Derive2(W, C) * 4.0;
	}
}
