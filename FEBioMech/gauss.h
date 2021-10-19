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

// gaussian quadrature
const int nint1 = 1;
const double gp1[nint1] = {
    0
};
const double gw1[nint1] = {
    2
};


const int nint2 = 2;
const double gp2[nint2] = {
    -sqrt(3.)/3.,
    sqrt(3.)/3.
};
const double gw2[nint2] = {
    1,
    1
};

const int nint3 = 3;
const double gp3[nint3] = {
    -sqrt(3./5.),
    0.,
    sqrt(3./5.)
};
const double gw3[nint3] = {
    5./9.,
    8./9.,
    5./9.
};

const int nint4 = 4;
const double gp4[nint4] = {
    -sqrt((3+2*sqrt(6./5.))/7.),
    -sqrt((3-2*sqrt(6./5.))/7.),
    sqrt((3-2*sqrt(6./5.))/7.),
    sqrt((3+2*sqrt(6./5.))/7.)
};
const double gw4[nint4] = {
    (18-sqrt(30.))/36.,
    (18+sqrt(30.))/36.,
    (18+sqrt(30.))/36.,
    (18-sqrt(30.))/36.
};

const int nint5 = 5;
const double gp5[nint5] = {
    -1./3.*sqrt(5+2*sqrt(10./7.)),
    -1./3.*sqrt(5-2*sqrt(10./7.)),
    0.,
    1./3.*sqrt(5-2*sqrt(10./7.)),
    1./3.*sqrt(5+2*sqrt(10./7.))
};
const double gw5[nint5] = {
    (322-13*sqrt(70.))/900.,
    (322+13*sqrt(70.))/900.,
    128./225.,
    (322+13*sqrt(70.))/900.,
    (322-13*sqrt(70.))/900.
};

const int nint6 = 6;
const double gp6[nint6] = {
    0.661209386466264,
    -0.661209386466264,
    -0.238619186083196,
    0.238619186083196,
    -0.932469514203152,
    0.932469514203152
};
const double gw6[nint6] = {
    0.360761573048138,
    0.360761573048138,
    0.467913934572691,
    0.467913934572691,
    0.171324492379170,
    0.171324492379170
};

const int nint7 = 7;
const double gp7[nint7] = {
    0.000000000000000,
    0.405845151377397,
    -0.405845151377397,
    -0.741531185599394,
    0.741531185599394,
    -0.949107912342758,
    0.949107912342758
};
const double gw7[nint7] = {
    0.417959183673469,
    0.381830050505118,
    0.381830050505118,
    0.279705391489276,
    0.279705391489276,
    0.129484966168869,
    0.129484966168869
};

const int nint8 = 8;
const double gp8[nint8] = {
    -0.183434642495649,
    0.183434642495649,
    -0.525532409916329,
    0.525532409916329,
    -0.796666477413626,
    0.796666477413626,
    -0.960289856497536,
    0.960289856497536
};
const double gw8[nint8] = {
    0.362683783378362,
    0.362683783378362,
    0.313706645877887,
    0.313706645877887,
    0.222381034453374,
    0.222381034453374,
    0.101228536290376,
    0.101228536290376
};

const int nint9 = 9;
const double gp9[nint9] = {
    0.000000000000000,
    -0.836031107326635,
    0.836031107326635,
    -0.968160239507626,
    0.968160239507626,
    -0.324253423403808,
    0.324253423403808,
    -0.613371432700590,
    0.613371432700590
};
const double gw9[nint9] = {
    0.330239355001259,
    0.180648160694857,
    0.180648160694857,
    0.081274388361574,
    0.081274388361574,
    0.312347077040002,
    0.312347077040002,
    0.260610696402935,
    0.260610696402935
};

const int nint10 = 10;
const double gp10[nint10] = {
    -0.148874338981631,
    0.148874338981631,
    -0.433395394129247,
    0.433395394129247,
    -0.679409568299024,
    0.679409568299024,
    -0.865063366688984,
    0.865063366688984,
    -0.973906528517171,
    0.973906528517171
};
const double gw10[nint10] = {
    0.295524224714752,
    0.295524224714752,
    0.269266719309996,
    0.269266719309996,
    0.219086362515982,
    0.219086362515982,
    0.149451349150580,
    0.149451349150580,
    0.0666713443086881,
    0.0666713443086881
};
