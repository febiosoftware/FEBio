#include "stdafx.h"
#include "FERemodelingElasticDomain.h"
#include "FERemodelingElasticMaterial.h"
#include <FECore/FEModel.h>
#include "FECore/FEAnalysis.h"

//-----------------------------------------------------------------------------
//! constructor
FERemodelingElasticDomain::FERemodelingElasticDomain(FEModel* pfem) : FEElasticSolidDomain(pfem)
{
}

//-----------------------------------------------------------------------------
// Reset data
void FERemodelingElasticDomain::Reset()
{
	FEElasticSolidDomain::Reset();

	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];
		int n = el.GaussPoints();
		for (int j=0; j<n; ++j) {
			FERemodelingMaterialPoint& pt = *el.GetMaterialPoint(j)->ExtractData<FERemodelingMaterialPoint>();
			pt.m_rhor = pme->Density();
		}
	}
}

//-----------------------------------------------------------------------------
//! \todo The material point initialization needs to move to the base class.
bool FERemodelingElasticDomain::Initialize(FEModel &fem)
{
	// initialize base class
	if (FEElasticSolidDomain::Initialize(fem) == false) return false;

	// get the elements material
	FEElasticMaterial* pme = m_pMat->GetElasticMaterial();
	for (size_t i=0; i<m_Elem.size(); ++i)
	{
		FESolidElement& el = m_Elem[i];

		// initialize referential solid density
		for (int n=0; n<el.GaussPoints(); ++n)
		{
			FERemodelingMaterialPoint& pt = *el.GetMaterialPoint(n)->ExtractData<FERemodelingMaterialPoint>();
			pt.m_rhor = pme->Density();
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
//! calculates the global stiffness matrix for this domain
void FERemodelingElasticDomain::StiffnessMatrix(FESolver* psolver)
{
	// repeat over all solid elements
	int NE = m_Elem.size();

	// I only need this for the element density stiffness
	FEModel& fem = psolver->GetFEModel();

	#pragma omp parallel for
	for (int iel=0; iel<NE; ++iel)
	{
		// element stiffness matrix
		matrix ke;
		vector<int> lm;
		
		FESolidElement& el = m_Elem[iel];

		// create the element's stiffness matrix
		int ndof = 3*el.Nodes();
		ke.resize(ndof, ndof);
		ke.zero();

		// calculate geometrical stiffness
		ElementGeometricalStiffness(el, ke);

		// calculate material stiffness
		ElementMaterialStiffness(el, ke);

		// calculate density stiffness
		ElementDensityStiffness(fem, el, ke);

		// assign symmetic parts
		// TODO: Can this be omitted by changing the Assemble routine so that it only
		// grabs elements from the upper diagonal matrix?
		for (int i=0; i<ndof; ++i)
			for (int j=i+1; j<ndof; ++j)
				ke[j][i] = ke[i][j];

		// get the element's LM vector
		UnpackLM(el, lm);

		// assemble element matrix in global stiffness matrix
		#pragma omp critical
		psolver->AssembleStiffness(el.m_node, lm, ke);
	}
}

//-----------------------------------------------------------------------------
//! calculates the solid element stiffness matrix
void FERemodelingElasticDomain::ElementStiffness(FEModel& fem, int iel, matrix& ke)
{
	FESolidElement& el = Element(iel);

	// calculate material stiffness (i.e. constitutive component)
	ElementMaterialStiffness(el, ke);

	// calculate geometrical stiffness
	ElementGeometricalStiffness(el, ke);

	// calculate density stiffness
	ElementDensityStiffness(fem, el, ke);

	// assign symmetic parts
	// TODO: Can this be omitted by changing the Assemble routine so that it only
	// grabs elements from the upper diagonal matrix?
	int ndof = 3*el.Nodes();
	int i, j;
	for (i=0; i<ndof; ++i)
		for (j=i+1; j<ndof; ++j)
			ke[j][i] = ke[i][j];
}


//-----------------------------------------------------------------------------
//! calculates element's density stiffness component for integration point n
//! \todo Remove the FEModel parameter. We only need to dt parameter, not the entire model.
//! \todo Problems seem to run better without this stiffness matrix
void FERemodelingElasticDomain::ElementDensityStiffness(FEModel& fem, FESolidElement &el, matrix &ke)
{
	int n, i, j;

	// make sure this is a remodeling material
	FERemodelingElasticMaterial* pmat = dynamic_cast<FERemodelingElasticMaterial*>(m_pMat); assert(pmat);
    
    const int NE = FEElement::MAX_NODES;
    vec3d gradN[NE];
    double *Grn, *Gsn, *Gtn;
    double Gr, Gs, Gt;
    vec3d kru, kur;
        
    // nr of nodes
    int neln = el.Nodes();
        
    // nr of integration points
    int nint = el.GaussPoints();
        
    // get the time step value
    double dt = fem.GetCurrentStep()->m_dt;
        
    // jacobian
    double Ji[3][3], detJt;
        
    // weights at gauss points
    const double *gw = el.GaussWeights();
        
    // density stiffness component for the stiffness matrix
    mat3d kab;
        
    // calculate geometrical element stiffness matrix
    for (n=0; n<nint; ++n)
    {
        // calculate jacobian
        double J = invjact(el, Ji, n);
        detJt = J*gw[n];
            
        Grn = el.Gr(n);
        Gsn = el.Gs(n);
        Gtn = el.Gt(n);
            
        for (i=0; i<neln; ++i)
        {
            Gr = Grn[i];
            Gs = Gsn[i];
            Gt = Gtn[i];
                
            // calculate global gradient of shape functions
            // note that we need the transposed of Ji, not Ji itself !
            gradN[i].x = Ji[0][0]*Gr+Ji[1][0]*Gs+Ji[2][0]*Gt;
            gradN[i].y = Ji[0][1]*Gr+Ji[1][1]*Gs+Ji[2][1]*Gt;
            gradN[i].z = Ji[0][2]*Gr+Ji[1][2]*Gs+Ji[2][2]*Gt;
        }
            
        // get the material point data
        FEMaterialPoint& mp = *el.GetMaterialPoint(n);
        double drhohat = pmat->m_pSupp->Tangent_Supply_Density(mp);
        mat3ds ruhat = pmat->m_pSupp->Tangent_Supply_Strain(mp);
        mat3ds crho = pmat->Tangent_Stress_Density(mp);
        double krr = (drhohat - 1./dt)/J;
            
        for (i=0; i<neln; ++i) {
            kur = (crho*gradN[i])/krr;
            for (j=0; j<neln; ++j)
            {
                kru = ruhat*gradN[j];
                kab = kur & kru;
                ke[3*i  ][3*j  ] -= kab(0,0)*detJt;
                ke[3*i  ][3*j+1] -= kab(0,1)*detJt;
                ke[3*i  ][3*j+2] -= kab(0,2)*detJt;
                    
                ke[3*i+1][3*j  ] -= kab(1,0)*detJt;
                ke[3*i+1][3*j+1] -= kab(1,1)*detJt;
                ke[3*i+1][3*j+2] -= kab(1,2)*detJt;
                    
                ke[3*i+2][3*j  ] -= kab(2,0)*detJt;
                ke[3*i+2][3*j+1] -= kab(2,1)*detJt;
                ke[3*i+2][3*j+2] -= kab(2,2)*detJt;
            }
        }
    }
}
