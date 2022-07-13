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
#include <FECore/FEBodyLoad.h>
#include <FECore/FEOctreeSearch.h>
#include <unordered_map>
#include "febiomix_api.h"

class FESolidElement;

class FEBIOMIX_API FESBMPointSource : public FEBodyLoad
{
public:
	FESBMPointSource(FEModel* fem);

	bool Init() override;

	void Update() override;

	void Accumulate(double dc);

	void SetPosition(const vec3d& pos);

	vec3d GetPosition() const;

	void SetSBMID(int id);

	int GetSBMID() const;

	void SetRate(double rate);

	void SetRadius(double radius);

	double GetRate() const;

	double GetdC() const;

	double GetdCp() const;

	void SetdC(double dC);

	void SetWeighVolume(bool b);

	void SetResetFlag(bool b);

	void SetAccumulateFlag(bool b);

	//std::vector<FEMaterialPoint*> FindIntInRadius();
	void FindIntInRadius(std::vector<FEMaterialPoint*> &possible_ints, double &total_elem);

	//! return all the elements in the given radius
	void FindNodesInRadius(std::vector<FEMaterialPoint*>& possible_ints, double& total_elem);

private:
	//void ResetSBM();

private:
	int		m_sbmId;	// The SBM ID that defins the cell's "concentration"
	vec3d	m_pos;	// the position (in reference coordinates)
	double	m_rate;	// density value at point source
	double	m_radius;
	double	m_Vc;
	bool	m_reset;
	bool	m_doReset;
	bool	m_weighVolume;
	bool	m_accumulate;	// accumulate species flag for the update
	double	m_dC = 0.0;		// total change of a species
	double	m_dCp = 0.0;

private:
	FEOctreeSearch		m_search;
	FESolidElement*		m_el;
	double				m_q[3];

	DECLARE_FECORE_CLASS();
};
