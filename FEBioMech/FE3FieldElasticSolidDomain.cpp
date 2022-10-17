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
void FE3FieldElasticSolidDomain::ELEM_DATA::Serialize(DumpStream& ar)
{
	ar & eJ;
	ar & ep;
	ar & Lk;
	ar & eJt;
	ar & eJp;
}

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
	m_Elem = d.m_Elem; 
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

	// allocate element data
	int NE = Elements();
	m_Data.resize(NE);

	// initialize element data
	for (int i=0; i<NE; ++i)
	{
		ELEM_DATA& d = m_Data[i];
		d.eJ = d.eJt = d.eJp = 1.0;
		d.ep = 0.0;
		d.Lk = 0.0;
	}

	return true;
}

//-----------------------------------------------------------------------------
void FE3FieldElasticSolidDomain::Reset()
{
	FEElasticSolidDomain::Reset();
	// initialize element data
	size_t NE = m_Data.size();
	for (size_t i=0; i<NE; ++i)
	{
		ELEM_DATA& d = m_Data[i];
        d.eJ = d.eJt = d.eJp = 1.0;
        d.ep = 0.0;
        d.Lk = 0.0;
	}
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FE3FieldElasticSolidDomain::PreSolveUpdate(const FETimeInfo& timeInfo)
{
    FEElasticSolidDomain::PreSolveUpdate(timeInfo);
    int NE = (int)m_Data.size();
#pragma omp parallel for
	for (int i=0; i<NE; ++i)
    {
        ELEM_DATA& d = m_Data[i];
        d.eJp = d.eJt;
    }
}

//-----------------------------------------------------------------------------
//! Stiffness matrix for three-field domain
void FE3FieldElasticSolidDomain::StiffnessMatrix(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();

	// repeat over all solid elements
	int NE = (int)m_Elem.size();
	#pragma omp parallel for
	for (int iel=0; iel<NE; ++iel)
	{
		FESolidElement& el = m_Elem[iel];

		// element stiffness matrix
		FEElementMatrix ke(el);

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate material stiffness (i.e. constitutive component)
		ElementMaterialStiffness(iel, ke);

		// calculate geometrical stiffness
		ElementGeometricalStiffness(iel, ke);

		// Calculate dilatational stiffness
		ElementDilatationalStiffness(fem, iel, ke);

		// assign symmetic parts
		// TODO: Can this be omitted by changing the Assemble routine so that it only
		// grabs elements from the upper diagonal matrix?
		for (int i=0; i<ndof; ++i)
			for (int j=i+1; j<ndof; ++j)
				ke[j][i] = ke[i][j];

		// get the element's LM vector
		vector<int> lm;
		UnpackLM(el, lm);
		ke.SetIndices(lm);

		// assemble element matrix in global stiffness matrix
		LS.Assemble(ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates dilatational element stiffness component for element iel

void FE3FieldElasticSolidDomain::ElementDilatationalStiffness(FEModel& fem, int iel, matrix& ke)
{
	int i, j, n;

	FESolidElement& elem = Element(iel);
	ELEM_DATA& ed = m_Data[iel];

	const int nint = elem.GaussPoints();
	const int neln = elem.Nodes();
	const int ndof = 3*neln;

	// get the elements material
	FEUncoupledMaterial* pmi = dynamic_cast<FEUncoupledMaterial*>(m_pMat);
	assert(pmi);

	// average global derivatives
	vector<double> gradN(3*neln);
	zero(gradN);

	// initial element volume
	double Ve = 0;

	// global derivatives of shape functions
	double Gx, Gy, Gz;
	const double *gw = elem.GaussWeights();

	// jacobian
	double Ji[3][3], Jt, J0;

	double *Gr, *Gs, *Gt;

	// repeat over gauss-points
	for (n=0; n<nint; ++n)
	{
		// calculate jacobian
		J0 = detJ0(elem, n);
		Jt = invjact(elem, Ji, n, m_alphaf)*m_alphaf;

		Jt *= gw[n];

		Ve += J0*gw[n];

		Gr = elem.Gr(n);
		Gs = elem.Gs(n);
		Gt = elem.Gt(n);

		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		for (i=0; i<neln; ++i)
		{
			Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
			Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
			Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

			gradN[3*i  ] += Gx*Jt;
			gradN[3*i+1] += Gy*Jt;
			gradN[3*i+2] += Gz*Jt;
		}
	}

	// get effective modulus
	double k = pmi->UJJ(ed.eJ);

	// next, we add the Lagrangian contribution
	// note that this term will always be zero if the material does not
	// use the augmented lagrangian
	k += ed.Lk*pmi->hpp(ed.eJ);

	// divide by initial volume
	k /= Ve;

	// calculate dilatational stiffness component
	// we only calculate the upper triangular part
	// since ke is symmetric.
	for (i=0; i<ndof; ++i)
		for (j=i; j<ndof; ++j)
			ke[i][j] += k*gradN[i]*gradN[j];
}

//-----------------------------------------------------------------------------
//! Calculates element material stiffness element matrix

void FE3FieldElasticSolidDomain::ElementMaterialStiffness(int iel, matrix &ke)
{
	FESolidElement &el = Element(iel);
	ELEM_DATA& ed = m_Data[iel];

	// Get the current element's data
	const int nint = el.GaussPoints();
	const int neln = el.Nodes();
	const int ndof = 3*neln;

	// global derivatives of shape functions
	const int NME = FEElement::MAX_NODES;
	vec3d G[NME];

	double Gxi, Gyi, Gzi;
	double Gxj, Gyj, Gzj;

	// The 'D' matrix
	double D[6][6] = {0};	// The 'D' matrix

	// The 'D*BL' matrix
	double DBL[6][3];

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// get the material
	FEUncoupledMaterial& mat = dynamic_cast<FEUncoupledMaterial&>(*m_pMat);

	// we need the following tensors for the dilational stiffness
	mat3dd I(1);
	tens4ds IxI = dyad1s(I);
	tens4ds I4  = dyad4s(I);
	tens4ds Cp = IxI - I4*2;

	// calculate element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian
		double detJt = ShapeGradient(el, n, G, m_alphaf)*gw[n]*m_alphaf;

		// setup the material point
		// NOTE: deformation gradient and determinant have already been evaluated in the stress routine
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// get the material's tangent
		// Note that we are only grabbing the deviatoric tangent. 
		// The other tangent terms depend on the pressure p
		// which we seperately
		tens4ds C = Cp*ed.ep + mat.DevTangent(mp);

		// get the 'D' matrix
		C.extract(D);

		// we only calculate the upper triangular part
		// since ke is symmetric. The other part is
		// determined below using this symmetry.
		for (int i = 0, i3 = 0; i<neln; ++i, i3 += 3)
		{
			Gxi = G[i].x;
			Gyi = G[i].y;
			Gzi = G[i].z;

			for (int j = i, j3 = i3; j<neln; ++j, j3 += 3)
			{
				Gxj = G[j].x;
				Gyj = G[j].y;
				Gzj = G[j].z;

				// calculate D*BL matrices
				DBL[0][0] = (D[0][0]*Gxj+D[0][3]*Gyj+D[0][5]*Gzj);
				DBL[0][1] = (D[0][1]*Gyj+D[0][3]*Gxj+D[0][4]*Gzj);
				DBL[0][2] = (D[0][2]*Gzj+D[0][4]*Gyj+D[0][5]*Gxj);

				DBL[1][0] = (D[1][0]*Gxj+D[1][3]*Gyj+D[1][5]*Gzj);
				DBL[1][1] = (D[1][1]*Gyj+D[1][3]*Gxj+D[1][4]*Gzj);
				DBL[1][2] = (D[1][2]*Gzj+D[1][4]*Gyj+D[1][5]*Gxj);

				DBL[2][0] = (D[2][0]*Gxj+D[2][3]*Gyj+D[2][5]*Gzj);
				DBL[2][1] = (D[2][1]*Gyj+D[2][3]*Gxj+D[2][4]*Gzj);
				DBL[2][2] = (D[2][2]*Gzj+D[2][4]*Gyj+D[2][5]*Gxj);

				DBL[3][0] = (D[3][0]*Gxj+D[3][3]*Gyj+D[3][5]*Gzj);
				DBL[3][1] = (D[3][1]*Gyj+D[3][3]*Gxj+D[3][4]*Gzj);
				DBL[3][2] = (D[3][2]*Gzj+D[3][4]*Gyj+D[3][5]*Gxj);

				DBL[4][0] = (D[4][0]*Gxj+D[4][3]*Gyj+D[4][5]*Gzj);
				DBL[4][1] = (D[4][1]*Gyj+D[4][3]*Gxj+D[4][4]*Gzj);
				DBL[4][2] = (D[4][2]*Gzj+D[4][4]*Gyj+D[4][5]*Gxj);

				DBL[5][0] = (D[5][0]*Gxj+D[5][3]*Gyj+D[5][5]*Gzj);
				DBL[5][1] = (D[5][1]*Gyj+D[5][3]*Gxj+D[5][4]*Gzj);
				DBL[5][2] = (D[5][2]*Gzj+D[5][4]*Gyj+D[5][5]*Gxj);

				ke[i3  ][j3  ] += (Gxi*DBL[0][0] + Gyi*DBL[3][0] + Gzi*DBL[5][0] )*detJt;
				ke[i3  ][j3+1] += (Gxi*DBL[0][1] + Gyi*DBL[3][1] + Gzi*DBL[5][1] )*detJt;
				ke[i3  ][j3+2] += (Gxi*DBL[0][2] + Gyi*DBL[3][2] + Gzi*DBL[5][2] )*detJt;

				ke[i3+1][j3  ] += (Gyi*DBL[1][0] + Gxi*DBL[3][0] + Gzi*DBL[4][0] )*detJt;
				ke[i3+1][j3+1] += (Gyi*DBL[1][1] + Gxi*DBL[3][1] + Gzi*DBL[4][1] )*detJt;
				ke[i3+1][j3+2] += (Gyi*DBL[1][2] + Gxi*DBL[3][2] + Gzi*DBL[4][2] )*detJt;

				ke[i3+2][j3  ] += (Gzi*DBL[2][0] + Gyi*DBL[4][0] + Gxi*DBL[5][0] )*detJt;
				ke[i3+2][j3+1] += (Gzi*DBL[2][1] + Gyi*DBL[4][1] + Gxi*DBL[5][1] )*detJt;
				ke[i3+2][j3+2] += (Gzi*DBL[2][2] + Gyi*DBL[4][2] + Gxi*DBL[5][2] )*detJt;
			}
		}
	}
}

//-----------------------------------------------------------------------------
//! calculates element's geometrical stiffness component for each integration point

void FE3FieldElasticSolidDomain::ElementGeometricalStiffness(int iel, matrix &ke)
{
	FESolidElement &el = Element(iel);

	const int NME = FEElement::MAX_NODES;
	double Gx[NME], Gy[NME], Gz[NME];
	double *Grn, *Gsn, *Gtn;
	double Gr, Gs, Gt;

	// nr of nodes
	int neln = el.Nodes();

	// nr of integration points
	int nint = el.GaussPoints();

	// jacobian
	double Ji[3][3], detJt;

	// weights at gauss points
	const double *gw = el.GaussWeights();

	// stiffness component for the initial stress component of stiffness matrix
	double kab;

	// calculate geometrical element stiffness matrix
	for (int n=0; n<nint; ++n)
	{
		// calculate jacobian
		detJt = invjact(el, Ji, n, m_alphaf)*gw[n]*m_alphaf;

		Grn = el.Gr(n);
		Gsn = el.Gs(n);
		Gtn = el.Gt(n);

		for (int i = 0; i<neln; ++i)
		{
			Gr = Grn[i];
			Gs = Gsn[i];
			Gt = Gtn[i];

			// calculate global gradient of shape functions
			// note that we need the transposed of Ji, not Ji itself !
			Gx[i] = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
			Gy[i] = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
			Gz[i] = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
		}

		// get the material point data
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());

		// element's Cauchy-stress tensor at gauss point n
		// s is the voight vector
		mat3ds& s = pt.m_s;

		for (int i = 0; i<neln; ++i)
			for (int j = i; j<neln; ++j)
			{
				kab = (Gx[i]*(s.xx()*Gx[j]+s.xy()*Gy[j]+s.xz()*Gz[j]) +
					   Gy[i]*(s.xy()*Gx[j]+s.yy()*Gy[j]+s.yz()*Gz[j]) + 
					   Gz[i]*(s.xz()*Gx[j]+s.yz()*Gy[j]+s.zz()*Gz[j]))*detJt;

				ke[3*i  ][3*j  ] += kab;
				ke[3*i+1][3*j+1] += kab;
				ke[3*i+2][3*j+2] += kab;
			}
	}
}

//-----------------------------------------------------------------------------
//! This function loops over all elements and updates the stress
void FE3FieldElasticSolidDomain::Update(const FETimeInfo& tp)
{
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

	// calculate the average dilatation and pressure
	double v = 0, vt = 0, V = 0;
	for (int n=0; n<nint; ++n)
	{
		v += detJt(el, n, m_alphaf)*gw[n];
        vt+= detJt(el, n)*gw[n];
		V += detJ0(el, n)*gw[n];
	}

	// calculate volume ratio
	ed.eJ = v / V;
    ed.eJt = vt / V;
    double eUt = mat.U(ed.eJt);
    double eUp = mat.U(ed.eJp);

	// Calculate pressure. This is a sum of a Lagrangian term and a penalty term
	//      <--- Lag. mult. -->  <-- penalty -->
	ed.ep = ed.Lk*mat.hp(ed.eJ) + mat.UJ(ed.eJ);
//	ed.ep = mat.UJ(ed.eJ);

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
        pt.m_p = ed.ep;

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

		L1 = L0 + k*pmi->h(ed.eJ);
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

			double hi = pmi->h(ed.eJ);
			ed.Lk += k*pmi->h(ed.eJ);
			ed.ep = ed.Lk*pmi->hp(ed.eJ) + k*log(ed.eJ)/ed.eJ;
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
