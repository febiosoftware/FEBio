#pragma once

#include "mat3d.h"

//-----------------------------------------------------------------------------
// The following classes are defined in this file
class tens3ds;	// symmetric 3o tensor
class tens3drs;	// right-conjugate symmetric 3o tensor
class tens3dls;	// left-conjugate symmetric 3o tensor
class tens3d;	// general 3o tensor (no symmetry)

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with full symmetry Tijk = Tjik = Tkji = Tikj = Tkij = Tjki (only 10 out of 27 components are unique)

// We store this tensor as a 1x10 array.
// [T] = [T111 T112 T113 T122 T123 T133 T222 T223 T233 T333]
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9

class tens3ds
{
public:
	enum { NNZ = 10 };

public:
	// constructors
	tens3ds();
	tens3ds(double m[10]);

	// arithmetic operators
	tens3ds operator + (const tens3ds& t) const;
	tens3ds operator - (const tens3ds& t) const;
	tens3ds operator * (double g) const;
	tens3ds operator / (double g) const;

	// arithmetic assignment operators
	tens3ds& operator += (const tens3ds& t);
	tens3ds& operator -= (const tens3ds& t);
	tens3ds& operator *= (double g);
	tens3ds& operator /= (double g);

	// unary operators
	tens3ds operator - () const;

	// initialize to zero
	void zero();

	vec3d contractdyad1(const vec3d& v);
	double tripledot(const tens3ds& H);
	
public:
	double d[NNZ];
};

tens3ds dyad3s(const vec3d& l, const vec3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with right-conjugate symmetry Gijk = Gikj (only 18 out of 27 components are unique)

// We store this tensor as a 1x10 array.
// [G] = [G111 G112 G113 G122 G123 G133 G211 G212 G213 G222 G223 G233 G311 G312 G313 G322 G323 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3drs
{
public:
	enum { NNZ = 18 };

public:
	// constructors
	tens3drs();
	tens3drs(double m[18]);

	// arithmetic operators
	tens3drs operator + (const tens3drs& t) const;
	tens3drs operator - (const tens3drs& t) const;
	tens3drs operator * (double g) const;
	tens3drs operator / (double g) const;

	// arithmetic assignment operators
	tens3drs& operator += (const tens3drs& t);
	tens3drs& operator -= (const tens3drs& t);
	tens3drs& operator *= (double g);
	tens3drs& operator /= (double g);

	// unary operators
	tens3drs operator - () const;

	// initialize to zero
	void zero();
	
	vec3d contractdyad1(const vec3d& v);
	vec3d contract2s(const mat3ds& s);
	double tripledot(const tens3drs& H);
	vec3d contractdyad2(const vec3d& v, const vec3d& w);
	tens3dls transpose();
	void contractleg2(const mat3d& F, int leg);

public:
	double d[NNZ];	// stored in column major order
};

tens3drs operator * (const mat3d& F, const tens3drs& t);
tens3drs dyad3rs(const vec3d& l, const vec3d& r);

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with left-conjugate symmetry Gijk = Gjik (only 18 out of 27 components are unique)

// We store this tensor as a 1x10 array.
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3dls
{
public:
	enum { NNZ = 18 };

public:
	// constructors
	tens3dls();
	tens3dls(double m[18]);

	// arithmetic operators
	tens3dls operator + (const tens3dls& t) const;
	tens3dls operator - (const tens3dls& t) const;
	tens3dls operator * (double g) const;
	tens3dls operator / (double g) const;

	tens3dls operator * (const mat3d& F) const;

	// arithmetic assignment operators
	tens3dls& operator += (const tens3dls& t);
	tens3dls& operator -= (const tens3dls& t);
	tens3dls& operator *= (double g);
	tens3dls& operator /= (double g);

	// unary operators
	tens3dls operator - () const;

	// initialize to zero
	void zero();

	// transpose
	tens3drs transpose();

public:
	double d[NNZ];	// stored in column major order
};

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with no symmetry (27 components)

// Due to symmetry we can store this tensor as a 1x27 array.
// [T] = [T111 T112 T113 T121 T122 T123 T131 T132 T133 T211 T212 T213 T221 T222 T223 T231 T232 T233 T311 T312 T313 T321 T322 T323 T331 T332 T333
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9  T10  T11  T12  T13  T14  T15  T16  T17  T18  T19  T20  T21  T22  T23  T24  T25  T26

class tens3d
{
public:
	enum { NNZ = 27 };

public:
	// constructors
	tens3d();
	tens3d(double m[10]);

	// arithmetic operators
	tens3d operator + (const tens3d& t) const;
	tens3d operator - (const tens3d& t) const;
	tens3d operator * (double g) const;
	tens3d operator / (double g) const;

	// arithmetic assignment operators
	tens3d& operator += (const tens3d& t);
	tens3d& operator -= (const tens3d& t);
	tens3d& operator *= (double g);
	tens3d& operator /= (double g);

	// unary operators
	tens3d operator - () const;

	// initialize to zero
	void zero();

	// return symmetric tens3ds
	tens3ds symm();
	
public:
	double d[NNZ];	// stored in column major order
};

tens3d operator + (const tens3dls& l, const tens3drs& r);
inline tens3d operator + (const tens3drs& r, const tens3dls& l) { return l+r; }

// The following file contains the actual definition of the class functions
#include "tens3ds.hpp"
#include "tens3drs.hpp"
#include "tens3dls.hpp"
#include "tens3d.hpp"
