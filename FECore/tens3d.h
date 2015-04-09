#pragma once

#include "mat3d.h"

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with full symmetry Tijk = Tjik = Tkji = Tikj = Tkij = Tjki (only 10 out of 27 components are unique)

// Due to symmetry we can store this tensor as a 1x10 array.
// [T] = [T111 T112 T113 T122 T123 T133 T222 T223 T233 T333]
//     =    T0   T1   T2   T3   T4   T5   T6   T7   T8   T9

class tens3ds
{
public:
	enum { NNZ = 10 };

	// default constructor
	tens3ds(){}
	
	tens3ds(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens3ds(double m[10])
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = m[i];
	}

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

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	void unit();

	vec3d contractdyad1(vec3d v);
	double tripledot3s(tens3ds H);
	
public:
	double d[NNZ];	// stored in column major order
};

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with right-conjugate symmetry Gijk = Gikj (only 18 out of 27 components are unique)

// Due to symmetry we can store this tensor as a 1x10 array.
// [G] = [G111 G112 G113 G122 G123 G133 G211 G212 G213 G222 G223 G233 G311 G312 G313 G322 G323 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3drs
{
public:
	enum { NNZ = 18 };

	// default constructor
	tens3drs(){}
	
	tens3drs(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens3drs(double m[18])
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = m[i];
	}

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

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	void unit();

	vec3d contractdyad1(vec3d v);
	vec3d contract2s(mat3ds s);
	vec3d contractdyad2(vec3d v, vec3d w);

	double tripledot3rs(tens3drs H);

public:
	double d[NNZ];	// stored in column major order
};

//-----------------------------------------------------------------------------
//! Class for 3rd order tensor with left-conjugate symmetry Gijk = Gjik (only 18 out of 27 components are unique)

// Due to symmetry we can store this tensor as a 1x10 array.
// [G] = [G111 G112 G113 G121 G122 G123 G131 G132 G133 G221 G222 G223 G231 G232 G233 G331 G332 G333]
//     =    G0   G1   G2   G3   G4   G5   G6   G7   G8   G9  G10  G11  G12  G13  G14  G15  G16  G17

class tens3dls
{
public:
	enum { NNZ = 18 };

	// default constructor
	tens3dls(){}
	
	tens3dls(const double g)
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = g;
	}

	tens3dls(double m[18])
	{
		for (int i = 0; i < NNZ; i++)
			d[i] = m[i];
	}

	// arithmetic operators
	tens3dls operator + (const tens3dls& t) const;
	tens3dls operator - (const tens3dls& t) const;
	tens3dls operator * (double g) const;
	tens3dls operator / (double g) const;

	// arithmetic assignment operators
	tens3dls& operator += (const tens3dls& t);
	tens3dls& operator -= (const tens3dls& t);
	tens3dls& operator *= (double g);
	tens3dls& operator /= (double g);

	// unary operators
	tens3dls operator - () const;

	// trace
	double tr() const;
	
	// initialize to zero
	void zero();

	void unit();

	tens3drs transpose();

	//vec3d contractdyad1(vec3d v);
	//vec3d contract2s(mat3ds s);
	
	//double tripledot3rs(tens3drs H);

public:
	double d[NNZ];	// stored in column major order
};

// The following file contains the actual definition of the class functions
#include "tens3ds.hpp"
#include "tens3drs.hpp"
#include "tens3dls.hpp"
