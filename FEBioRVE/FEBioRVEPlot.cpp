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
#include "stdafx.h"
#include <FECore/writeplot.h>
#include "FEBioRVEPlot.h"
#include "FEMicroMaterial.h"
#include "FEMicroMaterial2O.h"

//-----------------------------------------------------------------------------
//! Store the average deformation Hessian (G) for each element. 

class FEMicro2OG
{
public:
	tens3drs operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint2O& pt2O = *(mp.ExtractData<FEElasticMaterialPoint2O>());
		return pt2O.m_G;
	}
};

bool FEPlotElementGnorm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial2O* pme = dom.GetMaterial()->ExtractProperty<FEElasticMaterial2O>();
	if (pme == 0) return false;

	writeAverageElementValue<tens3drs, double>(dom, a, FEMicro2OG(), [](const tens3drs& m) { return m.tripledot(m); });

	return true;
}

//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress for each element.

class FEMicro1OPK1Stress
{
public:
	FEMicro1OPK1Stress(FEMicroMaterial* pm) : m_mat(pm) {}
	mat3d operator()(const FEMaterialPoint& mp)
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		FEMicroMaterialPoint* mmppt = mp_noconst.ExtractData<FEMicroMaterialPoint>();
		return m_mat->AveragedStressPK1(mmppt->m_rve, mp_noconst);
	}

private:
	FEMicroMaterial* m_mat;
};

class FEMicro2OPK1Stress
{
public:
	mat3d operator()(const FEMaterialPoint& mp)
	{
		FEMaterialPoint& mp_noconst = const_cast<FEMaterialPoint&>(mp);
		FEMicroMaterialPoint2O* mmppt = mp_noconst.ExtractData<FEMicroMaterialPoint2O>();
		return mmppt->m_rve.AveragedStressPK1(mp_noconst);
	}
};

bool FEPlotElementPK1norm::Save(FEDomain& dom, FEDataStream& a)
{
	FEMicroMaterial* pm1O = dynamic_cast<FEMicroMaterial*>(dom.GetMaterial());
	if (pm1O)
	{
		writeAverageElementValue<mat3d, double>(dom, a, FEMicro1OPK1Stress(pm1O), [](const mat3d& m) {return m.dotdot(m); });
		return true;
	}

	FEMicroMaterial2O* pm2O = dynamic_cast<FEMicroMaterial2O*>(dom.GetMaterial());
	if (pm2O == 0)
	{
		writeAverageElementValue<mat3d, double>(dom, a, FEMicro2OPK1Stress(), [](const mat3d& m) {return m.dotdot(m); });
		return true;
	}

	return false;
}


//-----------------------------------------------------------------------------
//! Store the norm of the average PK1 stress moment for each element. 

class FEMicro2OQK1
{
public:
	tens3drs operator()(const FEMaterialPoint& mp)
	{
		const FEElasticMaterialPoint2O& pt2O = *(mp.ExtractData<FEElasticMaterialPoint2O>());
		return pt2O.m_Q;
	}
};

bool FEPlotElementQK1norm::Save(FEDomain& dom, FEDataStream& a)
{
	FEElasticMaterial2O* pme = dynamic_cast<FEElasticMaterial2O*>(dom.GetMaterial());
	if (pme == 0) return false;

	// write solid element data
	writeAverageElementValue<tens3drs, double>(dom, a, FEMicro2OQK1(), [](const tens3drs& m) { return m.tripledot(m); });

	return true;
}

//-----------------------------------------------------------------------------
//! Element macro energy
bool FEPlotElementMicroEnergy::Save(FEDomain& dom, FEDataStream& a)
{
	FEMicroMaterial* pm1O = dynamic_cast<FEMicroMaterial*>(dom.GetMaterial());
	if (pm1O)
	{
		writeAverageElementValue<double>(dom, a, [](const FEMaterialPoint& mp) {
			const FEMicroMaterialPoint& mmpt = *(mp.ExtractData<FEMicroMaterialPoint>());
			return mmpt.m_micro_energy;
			});
		return true;
	}
	return false;
}
