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
#include "FEMeshDataInterpolator.h"

//! Helper class for mapping data between two point sets using moving least squares.
class FELeastSquaresInterpolator : public FEMeshDataInterpolator
{
	class Data
	{
	public:
		Data();
		Data(const Data& d);
		void operator = (const Data& d);

	public:
		matrix			A;
		std::vector<int>		index;
		std::vector<double>	W;
		std::vector<vec3d>	X;
		std::vector<int>		cpl;
	};

public:
	//! constructor
	FELeastSquaresInterpolator();

	//! Set dimension (2 or 3)
	void SetDimension(int d);

	//! Set the number of nearest neighbors to use (should be larger than 4)
	void SetNearestNeighborCount(int nnc);

	//! Set check for match flag. This will return the source value
	//! if the target point coincides
	void SetCheckForMatch(bool b);

	//! set the source points
	void SetSourcePoints(const std::vector<vec3d>& srcPoints);

	//! set the target points
	void SetTargetPoints(const std::vector<vec3d>& trgPoints);
	bool SetTargetPoint(const vec3d& trgPoint) override;

	//! initialize MLQ data
	bool Init() override;

	//! map source data onto target data
	//! input: sval - values of the source points
	//! output: tval - values at the target points
	bool Map(std::vector<double>& tval, std::function<double(int sourceNode)> src) override;

	// evaluate map
	double Map(int inode, std::function<double(int sourceNode)> src) override;
	vec3d MapVec3d(int inode, std::function<vec3d(int sourceNode)> src) override;

private:
	int		m_dim;
	int		m_nnc;
	bool	m_checkForMatch;
	std::vector<vec3d>	m_src;	// source points
	std::vector<vec3d>	m_trg;	// target points

	std::vector< Data >			m_data;
};
