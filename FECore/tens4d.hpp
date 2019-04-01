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

#pragma once
// NOTE: This file is automatically included from tens4d.h
// Users should not include this file manually!

inline double& tens4d::operator () (int i, int j, int k, int l)
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double tens4d::operator () (int i, int j, int k, int l) const
{
	const int m[3][3] = {{0,3,5},{6,1,4},{8,7,2}};
	const tens4d& T = (*this);
	return T(m[i][j], m[k][l]);
}

inline double& tens4d::operator () (int i, int j)
{
	const int m[9] = {0, 9, 18, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}

inline double tens4d::operator () (int i, int j) const
{
	const int m[9] = {0, 9, 18, 27, 36, 45, 54, 63, 72};
	return d[m[j]+i];
}
