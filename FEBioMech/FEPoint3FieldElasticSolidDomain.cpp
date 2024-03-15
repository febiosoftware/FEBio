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
#include "FE3FieldElasticSolidDomain.h"
#include "FEUncoupledMaterial.h"
#include <FECore/FEModel.h>
#include "FECore/log.h"

//-----------------------------------------------------------------------------
BEGIN_FECORE_CLASS(FE3FieldElasticSolidDomain, FEElasticSolidDomain)
	ADD_PARAMETER(m_blaugon, "laugon");
	ADD_PARAMETER(m_augtol , "atol");
	ADD_PARAMETER(m_naugmin, "minaug");
	ADD_PARAMETER(m_naugmax, "maxaug");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! constructor
FE3FieldElasticSolidDomain::FE3FieldElasticSolidDomain(FEModel* pfem) : FEElasticSolidDomain(pfem) 
{
	m_blaugon = false;
	m_augtol = 0.01;
	m_naugmin = 0;
	m_naugmax = 0;
}

//-----------------------------------------------------------------------------
//! \todo Do I really use this?
FE3FieldElasticSolidDomain& FE3FieldElasticSolidDomain::operator = (FE3FieldElasticSolidDomain& d) 
{ 
	m_pMesh = d.m_pMesh; 
	return (*this); 
}

//-----------------------------------------------------------------------------
bool FE3FieldElasticSolidDomain::DoAugmentations() const
{
	return m_blaugon;
}

//-----------------------------------------------------------------------------
//! Initialize the 3-field domain data
bool FE3FieldElasticSolidDomain::Init()
{
	// make sure the domain material uses an uncoupled formulation
	if (dynamic_cast<FEUncoupledMaterial*>(m_pMat) == 0) return false;
	if (FEElasticSolidDomain::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
//! This function loops over all elements and updates the stress
void FE3FieldElasticSolidDomain::Update(const FETimeInfo& tp)
{

	feLog("FE3FieldElasticSolidDomain Update ----- \n");
	bool berr = false;
	int NE = (int) m_Elem.size();
	#pragma omp parallel for shared(NE, berr)
	for (int i=0; i<NE; ++i)
	{
		try
		{
			UpdateElementStress(i, tp);
		}
		catch (NegativeJacobian e)
		{
			#pragma omp critical
			{
				berr = true;
				if (e.DoOutput()) feLogError(e.what());
			}
		}
	}

	if (berr) throw NegativeJacobianDetected();
}

//-----------------------------------------------------------------------------
//! This function updates the stresses for elements using the three-field formulation.
//! For such elements, the stress is a sum of a deviatoric stress, calculate by the
//! material and a dilatational term.
void FE3FieldElasticSolidDomain::UpdateElementStress(int iel, const FETimeInfo& tp)
{
    double dt = GetFEModel()->GetTime().timeIncrement;
    
	// get the material
	FEUncoupledMaterial& mat = *(dynamic_cast<FEUncoupledMaterial*>(m_pMat));

	// get the solid element
	FESolidElement& el = m_Elem[iel];
	ELEM_DATA& ed = m_Data[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// get the integration weights
	double* gw = el.GaussWeights();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	const int NME = FEElement::MAX_NODES;
	vec3d r0[NME], r[NME], vel[NME], acc[NME];
	for (int j=0; j<neln; ++j)
	{
        FENode& node = m_pMesh->Node(el.m_node[j]);
		r0[j] = node.m_r0;
        r[j] = node.m_rt*m_alphaf + node.m_rp*(1-m_alphaf);
        vel[j] = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2])*m_alphaf + node.m_vp*(1-m_alphaf);
        acc[j] = node.m_at*m_alpham + node.m_ap*(1-m_alpham);
	}


	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

        pt.m_p = pt.Lk*mat.hp(pt.J, pt.J_star) + mat.UJ(pt.J, pt.J_star);

		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are not used by most materials.
		mp.m_r0 = el.Evaluate(r0, n);
		mp.m_rt = el.Evaluate(r, n);

		// get the deformation gradient and determinant
        double Jt, Jp;
        mat3d Ft, Fp;
        Jt = defgrad(el, Ft, n);
        Jp = defgradp(el, Fp, n);
        pt.m_F = (m_alphaf==1.0? Ft : Ft*m_alphaf + Fp*(1-m_alphaf));
        pt.m_J = pt.m_F.det();
        mat3d Fi = pt.m_F.inverse();
        pt.m_L = (Ft - Fp)*Fi/dt;
        pt.m_v = el.Evaluate(vel, n);
        pt.m_a = el.Evaluate(acc, n);

        // update specialized material points
        m_pMat->UpdateSpecializedMaterialPoints(mp, tp);
        
		// calculate the stress at this material point
		// Note that we don't call the material's Stress member function.
		// The reason is that we need to use the averaged pressure for the element
		// and the Stress function uses the pointwise pressure. 
		// Therefore we call the DevStress function and add the pressure term
		// seperately. 
		pt.m_s = mat.DevStress(mp);
        
        // adjust stress for strain energy conservation
        if (m_alphaf == 0.5) 
		{
			// evaluate deviatoric strain energy at current and previous time
			mat3d Ftmp = pt.m_F;
			double Jtmp = pt.m_J;
			pt.m_F = Ft;
			pt.m_J = Jt;
			double Wt = mat.DevStrainEnergyDensity(mp);
			pt.m_F = Ftmp;
			pt.m_J = Jtmp;

			// store total strain energy density at current time
			pt.m_Wt = Wt + eUt;

			double Wp = pt.m_Wp;
            mat3ds D = pt.RateOfDeformation();
            double D2 = D.dotdot(D);
            if (D2 > 0)
                pt.m_s += D*(((Wt-Wp)/(dt*pt.m_J) - pt.m_s.dotdot(D))/D2);
            if (ed.eJt != ed.eJp)
                pt.m_s += mat3dd((eUt-eUp)/(ed.eJ*(ed.eJt-ed.eJp)));
        }
        else
            pt.m_s += mat3dd(ed.ep);
	}
}

//-----------------------------------------------------------------------------
//! Do augmentation
bool FE3FieldElasticSolidDomain::Augment(int naug)
{
	FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_pMat);
	assert(pmi);

	// make sure Augmented Lagrangian flag is on
	if (m_blaugon == false) return true;

	// do the augmentation
	int n;
	double normL0 = 0, normL1 = 0, L0, L1;
	double k = pmi->m_K;
	FEMesh& mesh = *m_pMesh;
	int NE = Elements();

	for (n=0; n<NE; ++n)
	{
		ELEM_DATA& ed = m_Data[n];

		L0 = ed.Lk;
		normL0 += L0*L0;

		L1 = L0 + k*pmi->h(ed.eJ, ed.eJ_star);
		normL1 += L1*L1;
	}

	normL0 = sqrt(normL0);
	normL1 = sqrt(normL1);

	// check convergence
	double pctn = 0;
	if (fabs(normL1) > 1e-10) pctn = fabs((normL1 - normL0)/normL1);

	feLog(" material %d\n", pmi->GetID());
	feLog("                        CURRENT         CHANGE        REQUIRED\n");
	feLog("   pressure norm : %15le%15le%15le\n", normL1, pctn, m_augtol);

	// check convergence
	bool bconv = true;
	if (pctn >= m_augtol) bconv = false;
	if (m_naugmin > naug) bconv = false;
	if ((m_naugmax > 0) && (m_naugmax <= naug)) bconv = true;

	// do the augmentation only if we have not yet converged
	if (bconv == false)
	{
		for (n=0; n<NE; ++n)
		{
			ELEM_DATA& ed = m_Data[n];

			double hi = pmi->h(ed.eJ, ed.eJ_star);
			ed.Lk += k*pmi->h(ed.eJ, ed.eJ_star);
			ed.ep = ed.Lk*pmi->hp(ed.eJ, ed.eJ_star) + k*log(ed.eJ/ed.eJ_star)/ed.eJ;
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//! Data serialization
void FE3FieldElasticSolidDomain::Serialize(DumpStream &ar)
{
	FEElasticSolidDomain::Serialize(ar);
	ar & m_Data;
}
