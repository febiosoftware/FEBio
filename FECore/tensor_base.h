#pragma once

//-----------------------------------------------------------------------------
// Template class for constructing some of the higher order tensor classes.
// Defines storage as well as some basic operations that do not depend on the 
// order in which the tensor components are stored.
template <class T, int NNZ> class tensor_base
{
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
template<class T, int NNZ> T tensor_base<T,NNZ>::operator + (const T& t) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] + t.d[i];
	return s;
}

// operator -
template<class T, int NNZ> T tensor_base<T,NNZ>::operator - (const T& t) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i] - t.d[i];
	return s;
}

// operator *
template<class T, int NNZ> T tensor_base<T,NNZ>::operator * (double g) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = g*d[i];
	return s;
}

// operator /
template<class T, int NNZ> T tensor_base<T,NNZ>::operator / (double g) const
{
	T s;
	for (int i=0; i<NNZ; i++) s.d[i] = d[i]/g;
	return s;
}

// assignment operator +=
template<class T, int NNZ> T& tensor_base<T,NNZ>::operator += (const T& t)
{
	for (int i=0; i<NNZ; i++) d[i] += t.d[i];
	return static_cast<T&>(*this);
}

// assignment operator -=
template<class T, int NNZ> T& tensor_base<T,NNZ>::operator -= (const T& t)
{
	for (int i=0; i<NNZ; i++) d[i] -= t.d[i];
	return static_cast<T&>(*this);
}

// assignment operator *=
template<class T, int NNZ> T& tensor_base<T,NNZ>::operator *= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] *= g;
	return static_cast<T&>(*this);
}

// assignment operator /=
template<class T, int NNZ> T& tensor_base<T,NNZ>::operator /= (double g)
{
	for (int i=0; i<NNZ; i++) d[i] /= g;
	return static_cast<T&>(*this);
}

// unary operator -
template<class T, int NNZ> T tensor_base<T,NNZ>::operator - () const
{
	T s;
	for (int i = 0; i < NNZ; i++) s.d[i] = -d[i];
	return s;
}

// intialize to zero
template<class T, int NNZ> void tensor_base<T,NNZ>::zero()
{
	for (int i = 0; i < NNZ; i++) d[i] = 0.0;
}
