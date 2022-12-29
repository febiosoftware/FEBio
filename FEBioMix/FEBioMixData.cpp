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
#include "FEBioMixData.h"
#include "FEBiphasicSolute.h"
#include "FETriphasic.h"
#include "FEMultiphasic.h"
#include "FECore/FEModel.h"
#include <FECore/FESolidDomain.h>

//-----------------------------------------------------------------------------
double FENodeConcentration::value(const FENode& node) 
{
	const int dof_C = GetFEModel()->GetDOFIndex("concentration", 0);
	return node.get(dof_C); 
}

//-----------------------------------------------------------------------------
double FELogElemFluidPressure::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_pa;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemFluidFluxZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FEBiphasicMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FEBiphasicMaterialPoint>();
		if (ppt) val += ppt->m_w.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteConcentration::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_ca[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[0].z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteRefConcentration::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_sbmr[0];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogFixedChargeDensity::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_cF;
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteConcentration_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_ca[m_nsol];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxX_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxY_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSoluteFluxZ_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_j[m_nsol].z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemElectricPotential::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_psi;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityX::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.x;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityY::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.y;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemCurrentDensityZ::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_Ie.z;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemSBMConcentration_::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i=0; i<nint; ++i)
	{
		FESolutesMaterialPoint* ppt = el.GetMaterialPoint(i)->ExtractData<FESolutesMaterialPoint>();
		if (ppt) val += ppt->m_sbmr[m_nsol];
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPorosity::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		const FEMaterialPoint* mp = el.GetMaterialPoint(i);
		const FEElasticMaterialPoint* et = (mp->ExtractData<FEElasticMaterialPoint>());
		const FEBiphasicMaterialPoint* pt = (mp->ExtractData<FEBiphasicMaterialPoint>());

		double p = (et && pt ? (1 - pt->m_phi0t / et->m_J) : 0.0);
		val += p;
	}
	return val / (double) nint;
}

//-----------------------------------------------------------------------------
double FELogElemPermeability::value(FEElement& el)
{
	FEDomain* dom = dynamic_cast<FEDomain*>(el.GetMeshPartition());
	if (dom == nullptr) return 0.0;
	FEMaterial* mat = dom->GetMaterial();
	if (mat == nullptr) return 0.0;

	FEBiphasic* biphasic = mat->ExtractProperty<FEBiphasic>();
	if (biphasic == nullptr) return 0.0;

	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		const FEMaterialPoint* mp = el.GetMaterialPoint(i);
		mat3ds K = biphasic->Permeability(const_cast<FEMaterialPoint&>(*mp));

		switch (m_comp)
		{
		case 0: val += K(0, 0); break;
		case 1: val += K(1, 1); break;
		case 2: val += K(2, 2); break;
		case 3: val += K(0, 1); break;
		case 4: val += K(1, 2); break;
		case 5: val += K(0, 2); break;
		}
	}
	return val / (double)nint;
}

//-----------------------------------------------------------------------------
double FELogElemSolidStress::value(FEElement& el)
{
	double val = 0.0;
	int nint = el.GaussPoints();
	for (int i = 0; i < nint; ++i)
	{
		const FEMaterialPoint* mp = el.GetMaterialPoint(i);
		const FEBiphasicMaterialPoint* pt = (mp->ExtractData<FEBiphasicMaterialPoint>());
		if (pt)
		{
			mat3ds ss = pt->m_ss;

			switch (m_comp)
			{
			case 0: val += ss(0, 0); break;
			case 1: val += ss(1, 1); break;
			case 2: val += ss(2, 2); break;
			case 3: val += ss(0, 1); break;
			case 4: val += ss(1, 2); break;
			case 5: val += ss(0, 2); break;
			}
		}
	}
	return val / (double)nint;
}

//=============================================================================

FELogDomainIntegralSBMConcentration::FELogDomainIntegralSBMConcentration(FEModel* fem, int sbm) : FELogDomainData(fem) 
{
	m_sbm = sbm;
}

double FELogDomainIntegralSBMConcentration::value(FEDomain& dom)
{
	double sum = 0.0;
	if (dynamic_cast<FESolidDomain*>(&dom))
	{
		FESolidDomain& solidDomain = dynamic_cast<FESolidDomain&>(dom);
		for (int i = 0; i < solidDomain.Elements(); ++i)
		{
			FESolidElement& el = solidDomain.Element(i);
			double val = 0.0;
			int nint = el.GaussPoints();
			double* gw = el.GaussWeights();
			for (int n = 0; n < nint; ++n)
			{
				FESolutesMaterialPoint* ppt = el.GetMaterialPoint(n)->ExtractData<FESolutesMaterialPoint>();
				if (ppt)
				{
					double Jw = solidDomain.detJt(el, n) * gw[n];
					val += ppt->m_sbmr[m_sbm] * Jw;
				}
			}

			sum += val;
		}
	}
	return sum;
}

//=============================================================================
FELogDomainIntegralSoluteConcentration::FELogDomainIntegralSoluteConcentration(FEModel* fem, int sol) : FELogDomainData(fem)
{
	m_nsol = sol;
}

double FELogDomainIntegralSoluteConcentration::value(FEDomain& dom)
{
	double sum = 0.0;
	if (dynamic_cast<FESolidDomain*>(&dom))
	{
		FESolidDomain& solidDomain = dynamic_cast<FESolidDomain&>(dom);
		for (int i = 0; i < solidDomain.Elements(); ++i)
		{
			FESolidElement& el = solidDomain.Element(i);
			double val = 0.0;
			int nint = el.GaussPoints();
			double* gw = el.GaussWeights();
			for (int n = 0; n < nint; ++n)
			{
				FESolutesMaterialPoint* ppt = el.GetMaterialPoint(n)->ExtractData<FESolutesMaterialPoint>();
				if (ppt)
				{
					double Jw = solidDomain.detJt(el, n) * gw[n];
					val += ppt->m_ca[m_nsol] * Jw;
				}
			}

			sum += val;
		}
	}
	return sum;
}
