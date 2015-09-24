#include "stdafx.h"
#include "FEElasticMultiscaleDomain1O.h"
#include "FEMicroMaterial.h"
#include "FECore/mat3d.h"
#include "FECore/tens6d.h"

//-----------------------------------------------------------------------------
//! constructor
FEElasticMultiscaleDomain1O::FEElasticMultiscaleDomain1O(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
}

//-----------------------------------------------------------------------------
//! Initialize element data
void FEElasticMultiscaleDomain1O::InitElements()
{
	const int NE = FEElement::MAX_NODES;
	vec3d x0[NE], xt[NE], r0, rt;
	FEMesh& m = *GetMesh();
	
	FEMicroMaterial* pmat = dynamic_cast<FEMicroMaterial*>(m_pMat);
	FEModel rve;
	rve.CopyFrom(pmat->m_mrve);
	rve.GetStep(0)->SetPrintLevel(FE_PRINT_NEVER);
	
	if (rve.Init() == false) throw FEMultiScaleException();

	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int neln = el.Nodes();
		for (int i=0; i<neln; ++i)
		{
			x0[i] = m.Node(el.m_node[i]).m_r0;
			xt[i] = m.Node(el.m_node[i]).m_rt;
		}

		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) 
		{
			r0 = el.Evaluate(x0, j);
			rt = el.Evaluate(xt, j);

			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
			FEMicroMaterialPoint& mmpt = *mp.ExtractData<FEMicroMaterialPoint>();
			pt.m_r0 = r0;
			pt.m_rt = rt;

			pt.m_J = defgrad(el, pt.m_F, j);
			mp.Init(false);

			if (mmpt.m_rve_init == false){
				mmpt.m_rve.CopyFrom(rve);
				mmpt.m_rve.Init();
				mmpt.m_rve_prev.CopyFrom(rve);
				mmpt.m_rve_prev.Init();
				mmpt.m_rve_init = true;}
		}
	}
}

//-----------------------------------------------------------------------------
//! Update element state data (mostly stresses, but some other stuff as well)
//! \todo Remove the remodeling solid stuff
void FEElasticMultiscaleDomain1O::UpdateElementStress(int iel, double dt)
{
	// get the solid element
	FESolidElement& el = m_Elem[iel];

	// get the number of integration points
	int nint = el.GaussPoints();

	// number of nodes
	int neln = el.Nodes();

	// nodal coordinates
	vec3d r0[FEElement::MAX_NODES];
	vec3d rt[FEElement::MAX_NODES];
	for (int j=0; j<neln; ++j)
	{
		r0[j] = m_pMesh->Node(el.m_node[j]).m_r0;
		rt[j] = m_pMesh->Node(el.m_node[j]).m_rt;
	}

	// get the integration weights
	double* gw = el.GaussWeights();

	// loop over the integration points and calculate
	// the stress at the integration point
	for (int n=0; n<nint; ++n)
	{
		FEMaterialPoint& mp = *el.GetMaterialPoint(n);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEMicroMaterialPoint& mmpt = *(mp.ExtractData<FEMicroMaterialPoint>());

		// material point coordinates
		// TODO: I'm not entirly happy with this solution
		//		 since the material point coordinates are used by most materials.
		pt.m_r0 = el.Evaluate(r0, n);
		pt.m_rt = el.Evaluate(rt, n);

		// get the deformation gradient and determinant
		pt.m_J = defgrad(el, pt.m_F, n);

		// calculate the stress at this material point
		FEMicroMaterial* pmat = dynamic_cast<FEMicroMaterial*>(m_pMat);
		
		// calculate the stress at this material point
		int plot_on = 0;
		int num_elem = Elements();

		// If it is a multi-element problem, plot for the middle integration point of each element 
		if (n == 13){
			plot_on = el.m_nID;}
		
		// If it is a multi-element problem, plot for the last integration point in the first and last element
		//if ((el.m_nID == 1 || el.m_nID == num_elem) && (n == nint-1)){
		//	plot_on = el.m_nID;}
		
		// If it is a single-element problem, plot for each intergration point of the element
		if (num_elem == 1){
			plot_on = num_elem;}

		pt.m_s = pmat->Stress1O(mp, plot_on, n+1);
	}
}