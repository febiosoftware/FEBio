/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
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
#include "tensor_base.h"

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class tens3ds;	// symmetric 3o tensor
class tens3drs;	// right-conjugate symmetric 3o tensor
class tens3dls;	// left-conjugate symmetric 3o tensor
class tens3d;	// general 3o tensor (no symmetry)

//-----------------------------------------------------------------------------
// traits for these classes defining the number of components
template <> class tensor_traits<tens3ds > {public: enum { NNZ = 10}; };
template <> class tensor_traits<tens3drs> {public: enum { NNZ = 18}; };
template <> class tensor_traits<tens3dls> {public: enum { NNZ = 18}; };
template <> class tensor_traits<tens3d  > {public: enum { NNZ = 27}; };

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with full symmetry Tijk = Tjik = Tkji = Tikj = Tkij = Tjki (only 10 out of 27 components are unique)

// We store this tensor as a 1x10 array.
// [T] = [T111 T112 T113 T122 T123 T133 T222 T223 T233 T333]
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9

class tens3ds : public tensor_base<tens3ds>
{
public:
	// constructors
	tens3ds(){}

	// access operator
	double operator () (int i, int j, int k) const;

	vec3d contractdyad1(const vec3d& v);
	double tripledot(const tens3ds& H);
};

tens3ds dyad3s(const vec3d& l, const vec3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with right-conjugate symmetry Gijk = Gikj (only 18 out of 27 components are unique)

// We store this tensor as a 1x18 array.
// [G] = [G111 G112 G113 G122 G123 G133 G211 G212 G213 G222 G223 G233 G311 G312 G313 G322 G323 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3drs : public tensor_base<tens3drs>
{
public:
	// constructors
	explicit tens3drs(double a);
	tens3drs(){}

	// access operator
	double  operator () (int i, int j, int k) const;
	double& operator () (int i, int j, int k);

	vec3d contractdyad1(const vec3d& v) const;
	vec3d contract2s(const mat3ds& s) const;
	double tripledot(const tens3drs& H) const;
	vec3d contractdyad2(const vec3d& v, const vec3d& w);
	tens3dls transpose();
	void contractleg2(const mat3d& F, int leg);
};

tens3drs operator * (const mat3d& F, const tens3drs& t);
tens3drs dyad3rs(const vec3d& l, const vec3d& r);
tens3drs dyad3rs(const mat3d& L, const vec3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with left-conjugate symmetry Gijk = Gjik (only 18 out of 27 components are unique)

// We store this tensor as a 1x18 array.
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3dls : public tensor_base<tens3dls>
{
public:
	// constructors
	tens3dls(){}

	tens3dls operator * (const mat3d& F) const;
    tens3dls operator * (const double& f) const;

	// transpose
	tens3drs transpose();
    vec3d trace();
    tens3d generalize();
};

tens3dls dyad3ls(const mat3ds& L, const vec3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with no symmetry (27 components)

// Due to symmetry we can store this tensor as a 1x27 array.
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26

class tens3d : public tensor_base<tens3d>
{
public:
	// constructors
	tens3d(){}
    explicit tens3d(double a);

	// access operators
	double operator () (int i, int j, int k) const;
	double& operator () (int i, int j, int k);

	// return symmetric tens3ds
	tens3ds symm();
    
    // right transpose
    tens3d transposer();
    
    //Contract by 2nd order tensor
    vec3d contract2(const mat3d& s) const;

    //Contract on right by vector
    mat3d contract1(const vec3d& v) const;
};

tens3d operator + (const tens3dls& l, const tens3drs& r);
inline tens3d operator + (const tens3drs& r, const tens3dls& l) { return l+r; }

// The following file contains the actual definition of the class functions
#include "tens3ds.hpp"
#include "tens3drs.hpp"
#include "tens3dls.hpp"
#include "tens3d.hpp"
