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

//-----------------------------------------------------------------------------
// traits class for tensors. Classes derived from tensor_base must specialize
// this class and define the NNZ enum variable which defines the number of components
// stored for that tensor class.
template <class T> class tensor_traits {};

//-----------------------------------------------------------------------------
// Template class for constructing some of the higher order tensor classes.
// Defines storage as well as some basic operations that do not depend on the 
// order in which the tensor components are stored.
template <class T> class tensor_base
{
	enum { NNZ = tensor_traits<T>::NNZ};

public:
	tensor_base(){}

	// arithmetic operators
	T operator + (const T& t) const;
	T operator - (const T& t) const;
	T operator * (double g) const;
	T operator / (double g) const;

	// arithmetic assignment operators
	T& operator += (const T& t);
	T& operator -= (const T& t);
	T& operator *= (double g);
	T& operator /= (double g);

	// unary operators
	T operator - () const;

	// initialize to zero
	void zero();

public:
	double	d[NNZ];
};

// operator +
template<class T> T tensor_base<T>::operator + (const T& t) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] + t.d[i];
	return s;
}

// operator -
template<class T> T tensor_base<T>::operator - (const T& t) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] - t.d[i];
	return s;
}

// operator *
template<class T> T tensor_base<T>::operator * (double g) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = g*d[i];
	return s;
}

// operator /
template<class T> T tensor_base<T>::operator / (double g) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i]/g;
	return s;
}

// assignment operator +=
template<class T> T& tensor_base<T>::operator += (const T& t)
{
	for (int i=0; i<NNZ; i++) d[i] += t.d[i];
	return static_cast<T&>(*this);
}

// assignment operator -=
template<class T> T& tensor_base<T>::operator -= (const T& t)
{
	for (int i=0; i<NNZ; i++) d[i] -= t.d[i];
	return static_cast<T&>(*this);
}

// assignment operator *=
template<class T> T& tensor_base<T>::operator *= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] *= g;
	return static_cast<T&>(*this);
}

// assignment operator /=
template<class T> T& tensor_base<T>::operator /= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] /= g;
	return static_cast<T&>(*this);
}

// unary operator -
template<class T> T tensor_base<T>::operator - () const
{
	T s;
	for (int i = 0; i < NNZ; i++) s.d[i] = -d[i];
	return s;
}

// intialize to zero
template<class T> void tensor_base<T>::zero()
{
	for (int i = 0; i < NNZ; i++) d[i] = 0.0;
}
