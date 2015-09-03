#pragma once
#include "FECore/FESolidDomain.h"
#include "FEBioMech/FEElasticDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D solid elements
//! Note that this class inherits from FEElasticSolidDomain since the biphasic domain
//! also needs to calculate elastic stiffness contributions.
//!
class FEBiphasicSolidDomain : public FESolidDomain, public FEElasticDomain
{
public:
	//! constructor
	FEBiphasicSolidDomain(FEModel* pfem);

	//! initialize class
	bool Initialize(FEModel& fem);

	//! activate
	void Activate();

	//! reset domain data
	void Reset();

	//! intitialize element data
	void InitElements();

	//! Unpack solid element data  (overridden from FEDomain)
	void UnpackLM(FEElement& el, vector<int>& lm);

	//! get the material (overridden from FEDomain)
	FEMaterial* GetMaterial() { return m_pMat; }

	//! set the material
	void SetMaterial(FEMaterial* pmat);

public:
	// update stresses (overridden from FEElasticDomain)
	void UpdateStresses(FEModel& fem);

	// update element stress
	void UpdateElementStress(int iel);

	//! calculates the global stiffness matrix for this domain
	void StiffnessMatrix(FESolver* psolver, bool bsymm, double dt);

	//! calculates the global stiffness matrix (steady-state case)
	void StiffnessMatrixSS(FESolver* psolver, bool bsymm, double dt);
	
public:
	// internal work (overridden from FEElasticDomain)
	void InternalForces(FEGlobalVector& R);

	//! internal fluid work
	void InternalFluidWork(vector<double>& R, double dt);

	//! internal fluid work (steady state analysis)
	void InternalFluidWorkSS(vector<double>& R, double dt);

public:
	//! element internal force vector
	void ElementInternalForce(FESolidElement& el, vector<double>& fe);
	
	//! Calculates the internal fluid forces
	bool ElementInternalFluidWork(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! Calculates the internal fluid forces for steady-state response
	bool ElementInternalFluidWorkSS(FESolidElement& elem, vector<double>& fe, double dt);
	
	//! calculates the element biphasic stiffness matrix
	bool ElementBiphasicStiffness(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the element biphasic stiffness matrix for steady-state response
	bool ElementBiphasicStiffnessSS(FESolidElement& el, matrix& ke, bool bsymm, double dt);
	
	//! calculates the solid element stiffness matrix
	void SolidElementStiffness(FESolidElement& el, matrix& ke);

	//! geometrical stiffness
	void ElementGeometricalStiffness(FESolidElement &el, matrix &ke);

	//! material stiffness component
	void ElementBiphasicMaterialStiffness(FESolidElement& el, matrix& ke);

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

protected:
	FEBiphasic*	m_pMat;
};
