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
#include <functional>
#include "fecore_api.h"

template <class T> T weightedAverage(T* d, double* w, int n)
{
	T s = d[0] * w[0];
	for (int i = 1; i < n; ++i) s += d[i] * w[i];
	return s;
}

template <class T> T weightedAverage(T* d, double* w, int n, std::function<T(const T&)> fnc)
{
	T s = fnc(d[0]) * w[0];
	for (int i = 1; i < n; ++i) s += fnc(d[i]) * w[i];
	return s;
}

FECORE_API mat3ds weightedAverageStructureTensor(mat3ds* d, double* w, int n);

// evaluate Log_p (X)
FECORE_API mat3ds Log(const mat3ds& p, const mat3ds& X);
FECORE_API mat3ds Exp(const mat3ds& p, const mat3ds& X);
