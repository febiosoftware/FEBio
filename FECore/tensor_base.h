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
