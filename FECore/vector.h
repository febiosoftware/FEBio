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
#include "vec3d.h"
#include "fecore_api.h"
using namespace std;

class FEMesh;

double FECORE_API operator*(const vector<double>& a, const vector<double>& b);
vector<double> FECORE_API operator - (vector<double>& a, vector<double>& b);
template<typename T> void zero(vector<T>& a) { fill(a.begin(), a.end(), T(0)); }
template<> inline void zero<vec3d>(vector<vec3d>& a) { fill(a.begin(), a.end(), vec3d(0,0,0)); }
template<typename T> void assign(vector<T>& a, const T& v) { fill(a.begin(), a.end(), v); }
void FECORE_API operator+=(vector<double>& a, const vector<double>& b);
void FECORE_API operator-=(vector<double>& a, const vector<double>& b);
void FECORE_API operator*=(vector<double>& a, double b);
vector<double> FECORE_API operator+(const vector<double>& a, const vector<double>& b);

// copy vector and scale
void FECORE_API vcopys(vector<double>& a, const vector<double>& b, double s);

// add scaled vector
void FECORE_API vadds(vector<double>& a, const vector<double>& b, double s);
void FECORE_API vsubs(vector<double>& a, const vector<double>& b, double s);

// vector subtraction: a = l - r
void FECORE_API vsub(vector<double>& a, const vector<double>& l, const vector<double>& r);

// scale each component of a vector
void FECORE_API vscale(vector<double>& a, const vector<double>& s);

// gather operation (copy mesh data to vector)
void FECORE_API gather(vector<double>& v, FEMesh& mesh, int ndof);
void FECORE_API gather(vector<double>& v, FEMesh& mesh, const vector<int>& dof);

// scatter operation (copy vector data to mesh)
void FECORE_API scatter(vector<double>& v, FEMesh& mesh, int ndof);
void FECORE_API scatter3(vector<double>& v, FEMesh& mesh, int ndof1, int ndof2, int ndof3);

// calculate l2 norm of vector
double FECORE_API l2_norm(const vector<double>& v);
double l2_norm(double* x, int n);

#endif // AFX_VECTOR_H__9F132D73_20B9_4AE9_A40B_EE4FB9D0FABD__INCLUDED_
