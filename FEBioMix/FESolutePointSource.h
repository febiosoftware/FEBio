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

class FEBIOMIX_API FESolutePointSource : public FEBodyLoad
{
public:
	FESolutePointSource(FEModel* fem);

	bool Init() override;

	void Accumulate(double dc);

	void Update() override;

	void SetPosition(const vec3d& v);

	vec3d GetPosition() const;

	void SetSoluteID(int soluteID);

	int GetSoluteID() const;

	void SetRate(double rate);

	void SetRadius(double radius);

	double GetRate() const;

	double GetdC() const;

	double GetdCp() const;

	void SetdC(double dC);

	void SetdCp(double dCp);

	void SetAccumulateFlag(bool b);

	void SetAccumulateCAFlag(bool b);

	//! Evaluate force vector
	void LoadVector(FEGlobalVector& R) override;

	//! evaluate stiffness matrix
	void StiffnessMatrix(FELinearSystem& S) override;

	//! return all the elements in the given radius
	void FindNodesInRadius(std::vector<FEElement*>& possible_nodes, double& total_elem);

	vec3d ClampNatC(double r[3]);

private:
	int		m_soluteId;	//!< solute ID
	double	m_rate;		//!< production rate
	vec3d	m_pos;		//!< position of source
	bool	m_accumulate = false; //!< accumulate flag
	bool	m_accumulate_ca; //! < accumulate actual concentration flag
	double	m_radius;
	double	m_Vc;
	double	m_dC = 0.0;
	double	m_dCp = 0.0;

private:
	FEOctreeSearch		m_search;
	FESolidElement*		m_el;
	double				m_q[3];
	int					m_dofC;

	DECLARE_FECORE_CLASS();
};
