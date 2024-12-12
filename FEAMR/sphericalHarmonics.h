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

#include <memory>
#include <FECore/matrix.h>
#include "feamr_api.h"

FEAMR_API void getSphereCoords(int numPts, const double* xCoords, const double* yCoords, const double* zCoords, double* theta, double* phi);

FEAMR_API std::unique_ptr<matrix> compSH(int order, int numPts, double* theta, double* phi);

FEAMR_API double harmonicY(int degree, int order, double theta, double phi, int numType);

FEAMR_API void reconstructODF(std::vector<double>& sphHarm, std::vector<double>& ODF, int numPts, double* theta, double* phi);

FEAMR_API void altGradient(int order, std::vector<double>& sphHarm, std::vector<double>& gradient);

FEAMR_API void remesh(std::vector<double>& gradient, double lengthScale, double hausd, double grad, std::vector<vec3d>& nodePos, std::vector<vec3i>& elems);
FEAMR_API void remeshFull(std::vector<double>& gradient, double lengthScale, double hausd, double grad, std::vector<vec3d>& nodePos, std::vector<vec3i>& elems);

// Taken from std::assoc_legendre definition in GCC
template<typename _Tp>
_Tp
__assoc_legendre_p(unsigned int __l, unsigned int __m, _Tp __x,
            _Tp __phase = _Tp(+1));
