/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in 
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
#include "vec3d.h"
#include "quatd.h"

//-----------------------------------------------------------------------------
// Class that defines an affine transformation (scale, rotate, translate).
// This currently applies the transformation as follows: 
// 1. scale : the scale is applied in the local coordinate system
// 2. rotate: rotation from local to global coordinates
// 3. translate: translate to a global position
// 
class FECORE_API FETransform
{
public:
	FETransform();

	// set the scale factors
	void SetScale(double sx, double sy, double sz);

	// set the translation
	void SetTranslation(const vec3d& t);

	// set the rotation quaternion
	void SetRotation(const quatd& q);

	// set the rotation vector (uses degrees)
	void SetRotation(const vec3d& r);

	// set rotation via Euler angles Tait-Bryan (Z,Y,X) convention (in degrees)
	void SetRotation(double X, double Y, double Z);

	// apply transformation
	vec3d Transform(const vec3d& r) const;

private:
	double	m_scl[3];	// scale factors
	vec3d	m_pos;		// translation (global space)
	quatd	m_rot;		// rotation
};
