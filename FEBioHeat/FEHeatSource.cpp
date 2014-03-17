#include "FEHeatSource.h"
#include "FECore/FEModel.h"

BEGIN_PARAMETER_LIST(FEHeatSource, FEBodyLoad);
	ADD_PARAMETER(m_Q, FE_PARAM_DOUBLE, "Q");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEHeatSource::FEHeatSource(FEModel* pfem) : FEBodyLoad(pfem)
{
	m_Q = 0;
}

//-----------------------------------------------------------------------------
void FEHeatSource::Residual(FEGlobalVector& R)
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEHeatSolidDomain* psd = dynamic_cast<FEHeatSolidDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			vector<double> fe;
			vector<int> lm;
			int NE = psd->Elements();
			for (int i=0; i<NE; ++i)
			{
				FESolidElement& el = psd->Element(i);
				int ne = el.Nodes();
				fe.resize(ne);
				ElementResidual(*psd, el, fe);
				psd->UnpackLM(el, lm);
				R.Assemble(el.m_node, lm, fe);
			}
		}
	}
}

//-----------------------------------------------------------------------------
void FEHeatSource::ElementResidual(FEHeatSolidDomain& dom, FESolidElement& el, vector<double>& fe)
{
	zero(fe);
	double* w = el.GaussWeights();
	int ne = el.Nodes();
	int ni = el.GaussPoints();
	for (int n=0; n<ni; ++n)
	{
		double* H = el.H(n);
		double J = dom.detJt(el, n);
		for (int i=0; i<ne; ++i)
		{
			fe[i] += m_Q*H[i]*J*w[n];
		}
	}
}

//-----------------------------------------------------------------------------
void FEHeatSource::Serialize(DumpFile& ar)
{
	FEBodyLoad::Serialize(ar);
}
