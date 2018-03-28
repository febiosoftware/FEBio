//
//  FEBiphasicShellDomain.hpp
//  FEBioMix
//
//  Created by Gerard Ateshian on 12/1/16.
//  Copyright © 2016 febio.org. All rights reserved.
//

#ifndef FEBiphasicShellDomain_hpp
#define FEBiphasicShellDomain_hpp

#include <FEBioMech/FESSIShellDomain.h>
#include "FEBiphasicDomain.h"
#include "FEBiphasic.h"

//-----------------------------------------------------------------------------
//! Domain class for biphasic 3D shell elements
//!
class FEBiphasicShellDomain : public FESSIShellDomain, public FEBiphasicDomain
{
public:
    //! constructor
    FEBiphasicShellDomain(FEModel* pfem);
    
    //! initialize class
	bool Init() override;
    
    //! activate
    void Activate();
    
    //! reset domain data
    void Reset();
    
    //! intitialize element data
    void PreSolveUpdate(const FETimeInfo& timeInfo);
    
    //! Unpack shell element data  (overridden from FEDomain)
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
    void ElementInternalForce(FEShellElement& el, vector<double>& fe);
    
    //! element internal force vector (stead-state case)
    void ElementInternalForceSS(FEShellElement& el, vector<double>& fe);
    
    //! calculates the element biphasic stiffness matrix
    bool ElementBiphasicStiffness(FEShellElement& el, matrix& ke, bool bsymm);
    
    //! calculates the element biphasic stiffness matrix for steady-state response
    bool ElementBiphasicStiffnessSS(FEShellElement& el, matrix& ke, bool bsymm);

public: // overridden from FEElasticDomain, but not all implemented in this domain
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf);
    void ElementBodyForce(FEBodyForce& BF, FEShellElement& el, vector<double>& fe);
    void InertialForces(FEGlobalVector& R, vector<double>& F) {}
    void StiffnessMatrix(FESolver* psolver) {}
    void BodyForceStiffness(FESolver* psolver, FEBodyForce& bf);
    void ElementBodyForceStiffness(FEBodyForce& BF, FEShellElement &el, matrix &ke);
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
    int					m_dofSX;
    int					m_dofSY;
    int					m_dofSZ;
};

#endif /* FEBiphasicShellDomain_hpp */
