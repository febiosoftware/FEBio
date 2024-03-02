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
#include "ad2.h"

double ad2::Evaluate(std::function<ad2::number(ad2::mat3ds& C)> W, const ::mat3ds& C)
{
	ad2::mat3ds dC(C);
	return W(dC).r;
}

::mat3ds ad2::Derive(std::function<ad2::number(ad2::mat3ds& C)> W, const ::mat3ds& C)
{
	ad2::mat3ds dC(C);
	double S[6] = { 0.0 };
	for (int i = 0; i < 6; ++i)
	{
		dC[i].d1 = 1;
		S[i] = W(dC).d1;
		dC[i].d1 = 0;
	}
	return ::mat3ds(S[0], S[2], S[5], 0.5 * S[1], 0.5 * S[4], 0.5 * S[3]);
}

::tens4ds ad2::Derive2(std::function<ad2::number(ad2::mat3ds& C)> W, const ::mat3ds& C)
{
	constexpr int l[6] = { 0, 2, 5, 1, 4, 3 };
	constexpr double w[6] = { 1.0, 1.0, 1.0, 0.5, 0.5, 0.5 };
	ad2::mat3ds dC(C);
	double D[6][6] = { 0 };
	for (int i = 0; i < 6; ++i)
		for (int j = 0; j < 6; ++j)
		{
			dC[l[i]].d1 = 1;
			dC[l[j]].d2 = 1;
			double ddW = W(dC).dd;
			dC[l[i]].d1 = 0;
			dC[l[j]].d2 = 0;

			D[i][j] = ddW*w[i]*w[j];
		}

	return ::tens4ds(D);
}
