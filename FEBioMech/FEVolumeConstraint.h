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
#include <FECore/FESurfaceConstraint.h>
#include <FECore/FESurface.h>

//-----------------------------------------------------------------------------
class FEVolumeSurface : public FESurface
{
public:
	//! constructor
	FEVolumeSurface(FEModel* fem);

	//! Initialization
	bool Init();

	//! copy data
	void CopyFrom(FEVolumeSurface& s);

	//! serialization
	void Serialize(DumpStream& ar);

public:
	double Volume();

public:
	double	m_Lp;	//!< Lagrange multipler pressure
	double	m_p;	//!< applied pressure (= Lp + eps*DV)
	double	m_V0;	//!< Initial volume
	double	m_Vt;	//!< current volume
};

//-----------------------------------------------------------------------------
// This class implements a constraint that tries to maintain the volume of the 
// enclosed space using an isochoric pressure.
class FEVolumeConstraint : public FESurfaceConstraint
{
public:
	//! constructor
	FEVolumeConstraint(FEModel* pfem);

	~FEVolumeConstraint();

	void Activate() override;
	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void Serialize(DumpStream& ar) override;
	void CopyFrom(FENLConstraint* plc) override;

	// update state
	void Reset() override;
	void Update() override;

	//! Unpack surface element data
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! build connectivity for matrix profile
	void BuildMatrixProfile(FEGlobalMatrix& M) override;

	// get the surface
	FESurface* GetSurface() override;

public:
	double EnclosedVolume() const;

	double Pressure() const;

public:
	double	m_eps;		//!< penalty parameter
	double	m_atol;		//!< augmented Lagrangian tolerance
	bool	m_blaugon;	//!< augmentation flag

private:
	bool	m_binit;	//!< flag indicating whether the constraint is initialized

	// degrees of freedom
	// (TODO: find a better way of defining this. 
	//        I don't want to have to do this in each class)
	int	m_dofX;
	int	m_dofY;
	int	m_dofZ;

	FEVolumeSurface* m_s;	//!< the bounding surface

	DECLARE_FECORE_CLASS();
};
