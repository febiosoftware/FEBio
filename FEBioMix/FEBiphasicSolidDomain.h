#pragma once
#include "FECore/FESolidDomain.h"
#include "FEBiphasicDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSolidDomain : public FESolidDomain, public FEBiphasicDomain
{
public:
	//! constructor
	FEBiphasicSolidDomain(FEModel* pfem);

	//! initialize class
	bool Initialize();

	//! activate
	void Activate();

	//! reset domain data
	void Reset();

	//! intitialize element data
	void PreSolveUpdate(const FETimeInfo& timeInfo);

	//! Unpack solid element data  (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

public:
	// update domain data
	void Update(const FETimeInfo& tp);

	// update element stress
	void UpdateElementStress(int iel);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm);

	//! calculates the global stiffness matrix (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm);
	
public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

    // internal work (steady-state case)
    void InternalForcesSS(FEGlobalVector& R);
    
public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	
    //! element internal force vector (steady-state case)
    void ElementInternalForceSS(FESolidElement& el, vector<double>& fe);
    
	//! calculates the element biphasic stiffness matrix
	bool ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm);
	
	//! calculates the element biphasic stiffness matrix for steady-state response
	bool ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm);
	
public: // overridden from FEElasticDomain, but not all implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf);
    void ElementBodyForce(FEBodyForce& BF, FESolidElement& el, vector<double>& fe);
	void InertialForces(FEGlobalVector& R, vector<double>& F) {}
	void StiffnessMatrix(FESolver* psolver) {}
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);
    void ElementBodyForceStiffness(FEBodyForce& BF, FESolidElement &el, matrix &ke);
	void MassMatrix(FESolver* psolver, double scale) {}

public: // biphasic domain "properties"
	// NOTE: I'm thinking about defining properties for domain classes. These would be similar to material
	// properties (and may require material properties to be evaluated), but are different in that they are
	// not meant to be customized. For example, the biphasic solver assumes Darcy's law in the evaluation 
	// of the fluid flux. Although one can see the fluid flux as a material property, since the code makes explicit
	// use of this constitutive equation (apparent from the fact that the biphasic material needs to define the permeability and its
	// strain derivate) it is not a true material property: i.e. it is not meant to be changed and is an inherent
	// assumption in this implementation. Consequently, the fluid flux would be a good example of a domain property.
	// That is why I've taken this calculation out of the FEBiphasic class and placed it here. 
	vec3d FluidFlux(FEMaterialPoint& mp);

	// Evaluate the nodal pressures
	// Note that the data vector stores the values for all of the nodes of the mesh, not just the domain nodes.
	// The values will be set to zero for nodes that don't belong to this domain.
	void GetNodalPressures(vector<double>& data);

private:
	// NOTE: This is a temporary construction. Just trying something out here.
	// Idea is here to construct a data export for the nodal pressures (the ones from the integration points, not the nodal dofs).
	vector<double>	m_nodePressure;	//!< nodal pressures projected from the integration points

	// This function updates the m_nodePressure variable
	void UpdateNodalPressures();
};
