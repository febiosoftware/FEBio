// vector.h: interface for the vector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_)
#define AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <memory.h>
#include <vector>
#include <algorithm>
using namespace std;

double operator*(const vector<double>& a, const vector<double>& b);
template<typename T> vector<T> operator - (vector<T>& a, vector<T>& b);
template<typename T> void zero(vector<T>& a) { fill(a.begin(), a.end(), T(0)); }
template<typename T> void assign(vector<T>& a, const T& v) { fill(a.begin(), a.end(), v); }
vector<double>& operator += (vector<double>& a, const vector<double>& b);
vector<double>& operator *= (vector<double>& a, double b);
vector<double> operator + (const vector<double>& a, const vector<double>& b);

// copy vector and scale
void vcopys(vector<double>& a, const vector<double>& b, double s);

#endif // AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_
