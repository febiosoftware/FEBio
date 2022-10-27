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
#include "FEJacobiSolidSolver.h"
#include "FERigidConnector.h"
#include "FESlidingElasticInterface.h"
#include "FE3FieldElasticSolidDomain.h"
#include "FE3FieldElasticShellDomain.h"
#include "FEElasticTrussDomain.h"
#include "FEBodyForce.h"
#include "FEResidualVector.h"
#include "FEUncoupledMaterial.h"
#include "FEContactInterface.h"
#include "FESSIShellDomain.h"
#include <FECore/log.h>
#include <FECore/DOFS.h>
#include <FECore/sys.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/vector.h>
#include "FESolidLinearSystem.h"
#include "FEBioMech.h"
#include "FESolidAnalysis.h"
#include "FETrussMaterial.h"
#include "FELinearTrussDomain.h"

//-----------------------------------------------------------------------------
// define the parameter list
BEGIN_FECORE_CLASS(FEJacobiSolidSolver, FESolver)
	BEGIN_PARAM_GROUP("Nonlinear solver");	// make sure this matches FENewtonSolver. 
		ADD_PARAMETER(m_Dtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "dtol"        );
		ADD_PARAMETER(m_Etol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "etol");
		ADD_PARAMETER(m_Rtol      , FE_RANGE_GREATER_OR_EQUAL(0.0), "rtol");
		ADD_PARAMETER(m_rhoi      , "rhoi"        );
		ADD_PARAMETER(m_alpha     , "alpha"       );
		ADD_PARAMETER(m_beta      , "beta"        );
		ADD_PARAMETER(m_gamma     , "gamma"       );
		ADD_PARAMETER(m_maxref    , FE_RANGE_GREATER_OR_EQUAL(0.0), "max_refs");
		ADD_PARAMETER(m_Rmin, FE_RANGE_GREATER_OR_EQUAL(0.0), "min_residual");
		ADD_PARAMETER(m_Rmax, FE_RANGE_GREATER_OR_EQUAL(0.0), "max_residual");
		ADD_PARAMETER(m_bdivreform, "diverge_reform");
		ADD_PARAMETER(m_breformAugment, "reform_augment");
		ADD_PARAMETER(m_bzero_diagonal, "check_zero_diagonal");
		ADD_PARAMETER(m_zero_tol, "zero_diagonal_tol");
		ADD_PARAMETER(m_breformtimestep, "reform_each_time_step");
		ADD_PARAMETER(m_bdoreforms          , "do_reforms"  );
	END_PARAM_GROUP();
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
//! FEJacobiSolidSolver Construction
//
FEJacobiSolidSolver::FEJacobiSolidSolver(FEModel* pfem) : FESolver(pfem), m_rigidSolver(pfem),\
m_dofU(pfem), m_dofV(pfem), m_dofSQ(pfem), m_dofRQ(pfem), m_dofSU(pfem), m_dofSV(pfem), m_dofSA(pfem)
{
	// default values
	m_Rtol = 0;	// deactivate residual convergence 
	m_Dtol = 0.001;
	m_Etol = 0.01;
	m_Rmin = 1.0e-20;
	m_Rmax = 0;	// not used if zero
	m_maxref = 15;
	m_bdivreform = true;
	m_breformAugment = true;
	m_bzero_diagonal = false;
	m_zero_tol = 0.0;
	m_breformtimestep = true;
	m_bforceReform = true;
	m_bdoreforms = false;

	m_niter = 0;
	m_nreq = 0;

	// default Newmark parameters (trapezoidal rule)
    m_rhoi = -2;
    m_alpha = m_alphaf = 1.0;
    m_alpham = 1.0;
	m_beta  = 0.25;
	m_gamma = 0.5;

	m_plinsolve = nullptr;
	m_pK = nullptr;

    // get the DOF indices
	// TODO: Can this be done in Init, since there is no error checking
	if (pfem)
	{
		m_dofU.AddVariable(FEBioMech::GetVariableName(FEBioMech::DISPLACEMENT));
		m_dofSQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ROTATION));
		m_dofRQ.AddVariable(FEBioMech::GetVariableName(FEBioMech::RIGID_ROTATION));
		m_dofV.AddVariable(FEBioMech::GetVariableName(FEBioMech::VELOCTIY));
		m_dofSU.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_DISPLACEMENT));
		m_dofSV.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_VELOCITY));
		m_dofSA.AddVariable(FEBioMech::GetVariableName(FEBioMech::SHELL_ACCELERATION));
	}
}

//-----------------------------------------------------------------------------
FEJacobiSolidSolver::~FEJacobiSolidSolver()
{

}

//-----------------------------------------------------------------------------
//! Return the rigid solver
FERigidSolver* FEJacobiSolidSolver::GetRigidSolver()
{
	return &m_rigidSolver;
}

//-----------------------------------------------------------------------------
//! Generate warnings if needed
void FEJacobiSolidSolver::SolverWarnings()
{
	FEModel& fem = *GetFEModel();

    // Generate warning if rigid connectors are used with symmetric stiffness
    if (m_msymm == REAL_SYMMETRIC) {
        for (int i=0; i<fem.NonlinearConstraints(); ++i)
        {
            FENLConstraint* plc = fem.NonlinearConstraint(i);
            FERigidConnector* prc = dynamic_cast<FERigidConnector*>(plc);
            if (prc) {
                feLogWarning("Rigid connectors require non-symmetric stiffness matrix.\nSet symmetric_stiffness flag to 0 in Control section.");
                break;
            }
        }
        
        // Generate warning if sliding-elastic contact is used with symmetric stiffness
		if (fem.SurfacePairConstraints() > 0)
        {
            // loop over all contact interfaces
			for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
            {
				FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
                FESlidingElasticInterface* pbw = dynamic_cast<FESlidingElasticInterface*>(pci);
                if (pbw) {
					feLogWarning("The sliding-elastic contact algorithm runs better with a non-symmetric stiffness matrix.\nYou may set symmetric_stiffness flag to 0 in Control section.");
                    break;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::CalculateMassMatrix()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	vector<double> dummy(m_Mi);
	FEGlobalVector Mi(fem, m_Mi, dummy);
	matrix me;
	vector <int> lm;
	vector <double> el_lumped_mass;

	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		// check whether it is a solid domain
		FEElasticSolidDomain* pbd = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
		if (pbd)  // it is an elastic solid domain
		{
			FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(pbd->GetMaterial());

			// loop over all the elements
			for (int iel = 0; iel < pbd->Elements(); ++iel)
			{
				FESolidElement& el = pbd->Element(iel);
				pbd->UnpackLM(el, lm);

				int nint = el.GaussPoints();
				int neln = el.Nodes();

				me.resize(neln, neln);
				me.zero();

				// create the element mass matrix
				for (int n = 0; n < nint; ++n)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(n);
					double d = pme->Density(mp);
					double detJ0 = pbd->detJ0(el, n) * el.GaussWeights()[n];

					double* H = el.H(n);
					for (int i = 0; i < neln; ++i)
						for (int j = 0; j < neln; ++j)
						{
							double kab = H[i] * H[j] * detJ0 * d;
							me[i][j] += kab;
						}
				}

				// reduce to a lumped mass vector and add up the total
				el_lumped_mass.assign(3 * neln, 0.0);
				for (int i = 0; i < neln; ++i)
				{
					for (int j = 0; j < neln; ++j)
					{
						double kab = me[i][j];
						el_lumped_mass[3 * i] += kab;
						el_lumped_mass[3 * i + 1] += kab;
						el_lumped_mass[3 * i + 2] += kab;
					}
				}

				// assemble element matrix into inv_mass vector 
				Mi.Assemble(el.m_node, lm, el_lumped_mass);
			} // loop over elements
		}
		else if (dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd)))
		{
			FEElasticShellDomain* psd = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(nd));
			FESolidMaterial* pme = dynamic_cast<FESolidMaterial*>(psd->GetMaterial());
			// loop over all the elements
			for (int iel = 0; iel < psd->Elements(); ++iel)
			{
				FEShellElement& el = psd->Element(iel);
				psd->UnpackLM(el, lm);

				// create the element's stiffness matrix
				FEElementMatrix ke(el);
				int neln = el.Nodes();
				int ndof = 6 * el.Nodes();
				ke.resize(ndof, ndof);
				ke.zero();

				// calculate inertial stiffness
				psd->ElementMassMatrix(el, ke, 1.0);

				// reduce to a lumped mass vector and add up the total
				el_lumped_mass.assign(ndof, 0.0);
				for (int i = 0; i < ndof; ++i)
				{
					for (int j = 0; j < ndof; ++j)
					{
						double kab = ke[i][j];
						el_lumped_mass[i] += kab;
					}
				}
				// assemble element matrix into inv_mass vector 
				Mi.Assemble(el.m_node, lm, el_lumped_mass);
			}
		}
		else if (dynamic_cast<FELinearTrussDomain*>(&mesh.Domain(nd)))
		{
			FELinearTrussDomain* ptd = dynamic_cast<FELinearTrussDomain*>(&mesh.Domain(nd));

			// loop over all the elements
			for (int iel = 0; iel < ptd->Elements(); ++iel)
			{
				FETrussElement& el = ptd->Element(iel);
				ptd->UnpackLM(el, lm);

				// create the element's stiffness matrix
				FEElementMatrix ke(el);
				int neln = el.Nodes();
				ke.resize(neln, neln);
				ke.zero();

				// calculate inertial stiffness
				ptd->ElementMassMatrix(el, ke);

				// reduce to a lumped mass vector and add up the total
				el_lumped_mass.assign(3 * neln, 0.0);
				for (int i = 0; i < neln; ++i)
				{
					for (int j = 0; j < neln; ++j)
					{
						double kab = ke[i][j];
						el_lumped_mass[3 * i] += kab;
						el_lumped_mass[3 * i + 1] += kab;
						el_lumped_mass[3 * i + 2] += kab;
					}
				}

				// assemble element matrix into inv_mass vector 
				Mi.Assemble(el.m_node, lm, el_lumped_mass);
			}
		}
		else if (dynamic_cast<FEElasticTrussDomain*>(&mesh.Domain(nd)))
		{
			FEElasticTrussDomain* ptd = dynamic_cast<FEElasticTrussDomain*>(&mesh.Domain(nd));

			// loop over all the elements
			for (int iel = 0; iel < ptd->Elements(); ++iel)
			{
				FETrussElement& el = ptd->Element(iel);
				ptd->UnpackLM(el, lm);

				// create the element's stiffness matrix
				FEElementMatrix ke(el);
				int neln = el.Nodes();
				ke.resize(neln, neln);
				ke.zero();

				// calculate inertial stiffness
				ptd->ElementMassMatrix(el, ke);

				// reduce to a lumped mass vector and add up the total
				el_lumped_mass.assign(3 * neln, 0.0);
				for (int i = 0; i < neln; ++i)
				{
					for (int j = 0; j < neln; ++j)
					{
						double kab = ke[i][j];
						el_lumped_mass[3 * i] += kab;
						el_lumped_mass[3 * i + 1] += kab;
						el_lumped_mass[3 * i + 2] += kab;
					}
				}

				// assemble element matrix into inv_mass vector 
				Mi.Assemble(el.m_node, lm, el_lumped_mass);
			}
		}
		else
		{
			//				return false;
		}
	}

	// we need the inverse of the lumped masses later
	// Also, make sure the lumped masses are positive.
	for (int i = 0; i < m_Mi.size(); ++i)
	{
		//		if (m_Mi[i] <= 0.0) return false;
		if (m_Mi[i] != 0.0) m_Mi[i] = 1.0 / m_Mi[i];
	}

	return true;
}

//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::AllocateLinearSystem()
{
	// Now that we have determined the equation numbers we can continue
	// with creating the stiffness matrix. First we select the linear solver
	// The stiffness matrix is created in CreateStiffness
	// Note that if a particular solver was requested in the input file
	// then the solver might already be allocated. That's way we need to check it.
	if (m_plinsolve == 0)
	{
		FEModel* fem = GetFEModel();
		FECoreKernel& fecore = FECoreKernel::GetInstance();
		m_plinsolve = fecore.CreateDefaultLinearSolver(fem);
		if (m_plinsolve == 0)
		{
			feLogError("Unknown solver type selected\n");
			return false;
		}
	}

	feLogInfo("Selecting linear solver %s", m_plinsolve->GetTypeStr());

	Matrix_Type mtype = MatrixType();
	SparseMatrix* pS = m_plinsolve->CreateSparseMatrix(mtype);
	if ((pS == 0) && (m_msymm == REAL_SYMMETRIC))
	{
		// oh, oh, something went wrong. It's probably because the user requested a symmetric matrix for a 
		// solver that wants a non-symmetric. If so, let's force a non-symmetric format.
		pS = m_plinsolve->CreateSparseMatrix(REAL_UNSYMMETRIC);

		if (pS)
		{
			// Problem solved! Let's inform the user.
			m_msymm = REAL_UNSYMMETRIC;
			feLogWarning("The matrix format was changed to non-symmetric since the selected linear solver does not support a symmetric format.");
		}
	}

	// if the sparse matrix is still zero, we have a problem
	if (pS == 0)
	{
		feLogError("The selected linear solver does not support the requested matrix format.\nPlease select a different linear solver.");
		return false;
	}

	// clean up the stiffness matrix if we have one
	if (m_pK) delete m_pK; m_pK = 0;

	// Create the stiffness matrix.
	// Note that this does not construct the stiffness matrix. This
	// is done later in the CreateStiffness routine.
	m_pK = new FEGlobalMatrix(pS);
	if (m_pK == 0)
	{
		feLogError("Failed allocating stiffness matrix.");
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Allocates and initializes the data structures used by the FEJacobiSolidSolver
//
bool FEJacobiSolidSolver::Init()
{
	// allocate data vectors
	m_R0.assign(m_neq, 0);
	m_R1.assign(m_neq, 0);
	m_ui.assign(m_neq, 0);
	m_Ui.assign(m_neq, 0);
	m_Ut.assign(m_neq, 0);
	m_Fd.assign(m_neq, 0);

	// allocate storage for the sparse matrix that will hold the stiffness matrix data
	// we let the linear solver allocate the correct type of matrix format
	if (AllocateLinearSystem() == false) return false;

	// Base class initialization and validation
	if (FESolver::Init() == false) return false;

	// set the create stiffness matrix flag
	m_breshape = true;

	FEModel& fem = *GetFEModel();

    if (m_rhoi == -1) {
        // Euler integration
        m_alpha = m_alphaf = m_alpham = 1.0;
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else if ((m_rhoi >= 0) && (m_rhoi <= 1)) {
        // Generalized-alpha integration (2nd order system)
        m_alpha = m_alphaf = 1.0/(1+m_rhoi);
        m_alpham = (2-m_rhoi)/(1+m_rhoi);
        m_beta = pow(1 + m_alpham - m_alphaf,2)/4;
        m_gamma = 0.5 + m_alpham - m_alphaf;
    }
    else {
        // for any other value of rhoi, use the user-defined alpha, beta, gamma parameters
        m_alphaf = m_alpham = m_alpha;
    }
    
	// allocate vectors
//	m_Fn.assign(m_neq, 0);
	m_Fr.assign(m_neq, 0);
	m_Ui.assign(m_neq, 0);
	m_Ut.assign(m_neq, 0);
	m_up.assign(m_neq, 0);
	m_Mi.assign(m_neq, 0.0);

	// we need to fill the total displacement vector m_Ut
	FEMesh& mesh = fem.GetMesh();
	gather(m_Ut, mesh, m_dofU[0]);
	gather(m_Ut, mesh, m_dofU[1]);
	gather(m_Ut, mesh, m_dofU[2]);
	gather(m_Ut, mesh, m_dofSQ[0]);
	gather(m_Ut, mesh, m_dofSQ[1]);
	gather(m_Ut, mesh, m_dofSQ[2]);
    gather(m_Ut, mesh, m_dofSU[0]);
    gather(m_Ut, mesh, m_dofSU[1]);
    gather(m_Ut, mesh, m_dofSU[2]);

    SolverWarnings();

	// calculate the inverse mass vector for the explicit analysis
	if (CalculateMassMatrix() == false)
	{
		feLogError("Failed building mass matrix.");
		return false;
	}

	// set the dynamic update flag only if we are running a dynamic analysis
	bool b = (fem.GetCurrentStep()->m_nanalysis == FESolidAnalysis::DYNAMIC ? true : false);
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEElasticSolidDomain* d = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
        FEElasticShellDomain* s = dynamic_cast<FEElasticShellDomain*>(&mesh.Domain(i));
		if (d) d->SetDynamicUpdateFlag(b);
        if (s) s->SetDynamicUpdateFlag(b);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Save data to dump file

void FEJacobiSolidSolver::Serialize(DumpStream& ar)
{
	FESolver::Serialize(ar);

	if (ar.IsShallow()) return;

	if (ar.IsSaving())
	{
		ar << m_neq;
		ar << m_maxref;
	}
	else
	{
		ar >> m_neq;
		ar >> m_maxref;

		// realloc data
		if (m_neq > 0)
		{
			m_R0.assign(m_neq, 0);
			m_R1.assign(m_neq, 0);
			m_ui.assign(m_neq, 0);
			m_Fd.assign(m_neq, 0);
		}
	}

	ar & m_nrhs;
	ar & m_niter;
	ar & m_nref & m_ntotref;
	ar & m_naug;
	ar & m_nreq;

	ar & m_alphaf;
	ar & m_alpham;

	ar & m_Ut & m_Ui;

	if (ar.IsLoading())
	{
//		m_Fn.assign(m_neq, 0);
		m_Fr.assign(m_neq, 0);
//		m_Ui.assign(m_neq, 0);
	}

	// serialize rigid solver
	m_rigidSolver.Serialize(ar);
}

//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::InitEquations()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FESolver::InitEquations() == false) return false;

	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq);
	if (neq == -1) return false; 
	else m_neq = neq;

	// Next, we add any Lagrange Multipliers
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* lmc = fem.NonlinearConstraint(i);
		if (lmc->IsActive())
		{
			m_neq += lmc->InitEquations(m_neq);
		}
	}
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
		if (spc->IsActive())
		{
			m_neq += spc->InitEquations(m_neq);
		}
	}

	// All initialization is done
	return true;
}


//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::InitEquations2()
{
	// First call the base class.
	// This will initialize all equation numbers, except the rigid body equation numbers
	if (FESolver::InitEquations2() == false) return false;

	// store the number of equations we currently have
	m_nreq = m_neq;

	// Next, we assign equation numbers to the rigid body degrees of freedom
	int neq = m_rigidSolver.InitEquations(m_neq);
	if (neq == -1) return false;
	else m_neq = neq;

	// Next, we add any Lagrange Multipliers
	FEModel& fem = *GetFEModel();
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* lmc = fem.NonlinearConstraint(i);
		if (lmc->IsActive())
		{
			m_neq += lmc->InitEquations(m_neq);
		}
	}
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
		if (spc->IsActive())
		{
			m_neq += spc->InitEquations(m_neq);
		}
	}

	// All initialization is done
	return true;
}

//-----------------------------------------------------------------------------
//! Update the kinematics of the model, such as nodal positions, velocities,
//! accelerations, etc.
void FEJacobiSolidSolver::UpdateKinematics(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// update rigid bodies
	m_rigidSolver.UpdateRigidBodies(m_Ui, ui);

	// total displacements
	vector<double> U(m_Ut.size());
	for (size_t i=0; i<m_Ut.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

	// update flexible nodes
	// translational dofs
	scatter(U, mesh, m_dofU[0]);
	scatter(U, mesh, m_dofU[1]);
	scatter(U, mesh, m_dofU[2]);
	// rotational dofs
	scatter(U, mesh, m_dofSQ[0]);
	scatter(U, mesh, m_dofSQ[1]);
	scatter(U, mesh, m_dofSQ[2]);
    // shell dofs
    scatter(U, mesh, m_dofSU[0]);
    scatter(U, mesh, m_dofSU[1]);
    scatter(U, mesh, m_dofSU[2]);

	// make sure the boundary conditions are fullfilled
	int nbcs = fem.BoundaryConditions();
	for (int i = 0; i<nbcs; ++i)
	{
		FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
		if (bc.IsActive()) bc.Update();
	}

	// enforce the linear constraints
	// TODO: do we really have to do this? Shouldn't the algorithm
	// already guarantee that the linear constraints are satisfied?
	FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
	if (LCM.LinearConstraints() > 0)
	{
		LCM.Update();
	}

	// Update the spatial nodal positions
	// Don't update rigid nodes since they are already updated
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        if (node.m_rid == -1) {
			node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
        }
        node.m_dt = node.m_d0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2])
        - node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
	}

	// update velocity and accelerations
	// for dynamic simulations
	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_nanalysis == FESolidAnalysis::DYNAMIC)
	{
		int N = mesh.Nodes();
		double dt = fem.GetTime().timeIncrement;
		double a = 1.0 / (m_beta*dt);
		double b = a / dt;
		double c = 1.0 - 0.5/m_beta;
		for (int i=0; i<N; ++i)
		{
			FENode& n = mesh.Node(i);
			n.m_at = (n.m_rt - n.m_rp)*b - n.m_vp*a + n.m_ap*c;
			vec3d vt = n.m_vp + (n.m_ap*(1.0 - m_gamma) + n.m_at*m_gamma)*dt;
			n.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vt);
            
            // shell kinematics
            vec3d qt = n.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
            vec3d qp = n.get_vec3d_prev(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
            vec3d vqp = n.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
            vec3d aqp = n.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
            vec3d aqt = (qt - qp)*b - vqp*a + aqp*c;
            vec3d vqt = vqp + (aqp*(1.0 - m_gamma) + aqt*m_gamma)*dt;
            n.set_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2], aqt);
            n.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vqt);
        }
    }

	// update nonlinear constraints (needed for updating Lagrange Multiplier)
	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* nlc = fem.NonlinearConstraint(i);
		if (nlc->IsActive()) nlc->Update(m_Ui, ui);
	}
	for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
	{
		FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
		if (spc->IsActive()) spc->Update(ui);
	}
}

//-----------------------------------------------------------------------------
//! Update DOF increments
void FEJacobiSolidSolver::UpdateIncrements(vector<double>& Ui, vector<double>& ui, bool emap)
{
	FEModel& fem = *GetFEModel();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();
    
	// update rigid bodies
	m_rigidSolver.UpdateIncrements(Ui, ui, emap);
        
	// update flexible nodes
	int n;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
        
		// displacement dofs
		// current position = initial + total at prev conv step + total increment so far + current increment
		if ((n = node.m_ID[m_dofU[0]]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[m_dofU[1]]) >= 0) Ui[n] += ui[n];
		if ((n = node.m_ID[m_dofU[2]]) >= 0) Ui[n] += ui[n];
        
        // rotational dofs
        if ((n = node.m_ID[m_dofSQ[0]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSQ[1]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSQ[2]]) >= 0) Ui[n] += ui[n];
        
        // shell dofs
        if ((n = node.m_ID[m_dofSU[0]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSU[1]]) >= 0) Ui[n] += ui[n];
        if ((n = node.m_ID[m_dofSU[2]]) >= 0) Ui[n] += ui[n];
	}

	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc && plc->IsActive()) plc->UpdateIncrements(Ui, ui);
	}
}

//-----------------------------------------------------------------------------
//! Updates the current state of the model
void FEJacobiSolidSolver::Update(vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    tp.currentIteration = m_niter;
    
    // update EAS
    UpdateEAS(ui);
    UpdateIncrementsEAS(ui, true);

	// update kinematics
	UpdateKinematics(ui);

	// update model state
	fem.Update();
}

//-----------------------------------------------------------------------------
//! Update EAS
void FEJacobiSolidSolver::UpdateEAS(vector<double>& ui)
{
	FEModel& fem = *GetFEModel();

    FEMesh& mesh = fem.GetMesh();

    // update EAS on shell domains
    for (int i=0; i<mesh.Domains(); ++i) {
        FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
        if (sdom && sdom->IsActive()) sdom->UpdateEAS(ui);
    }
}

//-----------------------------------------------------------------------------
//! Update EAS
void FEJacobiSolidSolver::UpdateIncrementsEAS(vector<double>& ui, const bool binc)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
    
    // update EAS on shell domains
    for (int i=0; i<mesh.Domains(); ++i) {
        FESSIShellDomain* sdom = dynamic_cast<FESSIShellDomain*>(&mesh.Domain(i));
        if (sdom && sdom->IsActive()) sdom->UpdateIncrementsEAS(ui, binc);
    }
}

//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::InitStep(double time)
{
	FEModel& fem = *GetFEModel();

	// set time integration parameters
	FETimeInfo& tp = fem.GetTime();
	tp.alpha = m_alpha;
	tp.beta = m_beta;
	tp.gamma = m_gamma;
	tp.alphaf = m_alphaf;
	tp.alpham = m_alpham;

    // evaluate load curve values at current (or intermediate) time
	double t = tp.currentTime;
//	double dt = tp.timeIncrement;
//	double ta = (t > 0) ? t - (1-m_alpha)*dt : m_alpha*dt;
//	return FESolver::InitStep(ta);
    return FESolver::InitStep(t);
}

//-----------------------------------------------------------------------------
//! Prepares the data for the first BFGS-iteration. 
bool FEJacobiSolidSolver::PrepStep()
{
	FEModel& fem = *GetFEModel();

    FETimeInfo& tp = fem.GetTime();
	double dt = tp.timeIncrement;
	tp.augmentation = 0;
    
	// zero total displacements
	zero(m_Ui);

	// store previous mesh state
	// we need them for velocity and acceleration calculations
	FEMesh& mesh = fem.GetMesh();
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& ni = mesh.Node(i);
		ni.m_rp = ni.m_rt;
		ni.m_vp = ni.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
		ni.m_ap = ni.m_at;
        ni.m_dp = ni.m_dt;
		ni.UpdateValues();

        // initial guess at start of new time step
        // solid
        ni.m_at = ni.m_ap*(1-0.5/m_beta) - ni.m_vp/(m_beta*dt);
        vec3d vs = ni.m_vp + (ni.m_at*m_gamma + ni.m_ap*(1-m_gamma))*dt;
        ni.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], vs);
        
        // solid shell
        vec3d aqp = ni.get_vec3d_prev(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
        vec3d vqp = ni.get_vec3d_prev(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
        vec3d aqt = aqp*(1-0.5/m_beta) - vqp/(m_beta*dt);
        ni.set_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2], aqt);
        vec3d vqt = vqp + (aqt*m_gamma + aqp*(1-m_gamma))*dt;
        ni.set_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2], vqt);
    }

    // apply concentrated nodal forces
	// since these forces do not depend on the geometry
	// we can do this once outside the NR loop.
//	vector<double> dummy(m_neq, 0.0);
//	zero(m_Fn);
//	FEResidualVector Fn(*GetFEModel(), m_Fn, dummy);
//	NodalLoads(Fn, tp);

	// apply boundary conditions
	// we save the prescribed displacements increments in the ui vector
	vector<double>& ui = m_ui;
	zero(ui);
	zero(m_up);
	int nbc = fem.BoundaryConditions();
	for (int i=0; i<nbc; ++i)
	{
		FEBoundaryCondition& dc = *fem.BoundaryCondition(i);
		if (dc.IsActive()) dc.PrepStep(ui);
	}

	// do the linear constraints
	fem.GetLinearConstraintManager().PrepStep();

	// initialize rigid bodies
	m_rigidSolver.PrepStep(tp, ui);

	// intialize material point data
	// NOTE: do this before the stresses are updated
	// TODO: does it matter if the stresses are updated before
	//       the material point data is initialized
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		FEDomain& dom = mesh.Domain(i);
		if (dom.IsActive()) dom.PreSolveUpdate(tp);
	}

	// update model state
	fem.Update();

	for (int i = 0; i < fem.NonlinearConstraints(); ++i)
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc && plc->IsActive()) plc->PrepStep();
	}

	// see if we need to do contact augmentations
	m_baugment = false;
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface& ci = dynamic_cast<FEContactInterface&>(*fem.SurfacePairConstraint(i));
		if (ci.IsActive() && (ci.m_laugon == 1)) m_baugment = true;
	}

	// see if we need to do incompressible augmentations
	// TODO: Should I do these augmentations in a nlconstraint class instead?
	int ndom = mesh.Domains();
	for (int i = 0; i < ndom; ++i)
	{
		FEDomain* dom = &mesh.Domain(i);
		FE3FieldElasticSolidDomain* dom3f = dynamic_cast<FE3FieldElasticSolidDomain*>(dom);
		if (dom3f && dom3f->DoAugmentations()) m_baugment = true;

		FE3FieldElasticShellDomain* dom3fs = dynamic_cast<FE3FieldElasticShellDomain*>(dom);
		if (dom3fs && dom3fs->DoAugmentations()) m_baugment = true;
	}

	// see if we have to do nonlinear constraint augmentations
	if (fem.NonlinearConstraints() != 0) m_baugment = true;

	// see if we reform at the start of every time step
	bool breform = m_breformtimestep;

	// if the force reform flag was set, we force a reform
	// (This will be the case for the first time this is called, or when the previous time step failed)
	if (m_bforceReform)
	{
		breform = true;

		m_bforceReform = false;
	}

	// do the reform
	// NOTE: It is important for JFNK that the matrix is reformed before the 
	//       residual is evaluated, so do not switch these two calculations!
	if (breform)
	{
		// do the first stiffness formation
		if (ReformStiffness() == false) return false;
	}

	// calculate initial residual
	if (Residual(m_R0) == false) return false;

	// add the contribution from prescribed dofs
	m_R0 += m_Fd;

	// TODO: I can check here if the residual is zero.
	// If it is than there is probably no force acting on the system
	// if (m_R0*m_R0 < eps) bconv = true;

	//	double r0 = m_R0*m_R0;

	// do callback (we do it here since we want the RHS to be formed as well)
	GetFEModel()->DoCallback(CB_MATRIX_REFORM);

	return true;
}

//-----------------------------------------------------------------------------
void FEJacobiSolidSolver::SolveEquations(std::vector<double>& u, std::vector<double>& R)
{
	int neq = m_neq;
	SparseMatrix* K = *m_pK; // stiffness matrix

	// multiply K with u_p
	vector<double> Ku(neq);
	K->mult_vector(m_up.data(), Ku.data());

	// evaluate new update
	const FETimeInfo& ti = GetFEModel()->GetTime();
	double dt = ti.timeIncrement;

	double s = m_beta * dt * dt / m_alpham;
	for (int i = 0; i < neq; ++i) u[i] = s * m_Mi[i] * (R[i] - Ku[i]);

	// store solution for next iteration
	m_up = u;
}

//-----------------------------------------------------------------------------
// Performs the quasi-newton iterations.
bool FEJacobiSolidSolver::SolveStep()
{
	// initialize counters
	m_niter = 0;	// nr of iterations
	m_nrhs = 0;		// nr of RHS evaluations
	m_nref = 0;		// nr of stiffness reformations
	m_ntotref = 0;
	m_naug = 0;		// nr of augmentations

	vector<double> u0(m_neq);
	vector<double> Rold(m_neq);

	// convergence norms
	double	normR1;		// residual norm
	double	normE1;		// energy norm
	double	normU;		// displacement norm
	double	normu;		// displacement increment norm
	double	normRi;		// initial residual norm
	double	normEi;		// initial energy norm
	double	normEm;		// max energy norm
	double	normUi;		// initial displacement norm

	FEModel& fem = *GetFEModel();

	// Get the current step
	FEAnalysis* pstep = fem.GetCurrentStep();

	// set the time information
	FETimeInfo& tp = fem.GetTime();

	// prepare for the first iteration
	if (PrepStep() == false) return false;

	// loop until converged or when max nr of reformations reached
	bool bconv = false;		// convergence flag
	do
	{
		feLog(" %d\n", m_niter+1);

		// assume we'll converge. 
		bconv = true;

		// solve the equations
		SolveEquations(m_ui, m_R0);

		// do the line search
		// Update geometry
		Update(m_ui);

		// calculate residual at this point
		Residual(m_R1);
		double s = 1;

		// set initial convergence norms
		if (m_niter == 0)
		{
			normRi = fabs(m_R0*m_R0);
			normEi = fabs(m_ui*m_R0);
			normUi = fabs(m_ui*m_ui);
			normEm = normEi;
		}

		// calculate actual displacement increment
		// NOTE: We don't apply the line search directly to m_ui since we need the unscaled search direction for the QN update below
		int neq = (int)m_Ui.size();
		vector<double> ui(m_ui);
		for (int i = 0; i<neq; ++i) ui[i] *= s;

		// update total displacements
		UpdateIncrements(m_Ui, ui, false);

		// calculate norms
		normR1 = m_R1*m_R1;
		normu  = ui*ui;
		normU  = m_Ui*m_Ui;
		normE1 = fabs(ui*m_R1);

		// check for nans
		if (ISNAN(normR1) || ISNAN(normu)) throw NANDetected();

		// check residual norm
		if ((m_Rtol > 0) && (normR1 > m_Rtol*normRi)) bconv = false;	

		// check displacement norm
		if ((m_Dtol > 0) && (normu  > (m_Dtol*m_Dtol)*normU )) bconv = false;

		// check energy norm
		if ((m_Etol > 0) && (normE1 > m_Etol*normEi)) bconv = false;

		// check linestep size
//		if ((m_lineSearch->m_LStol > 0) && (s < m_lineSearch->m_LSmin)) bconv = false;

		// check energy divergence
		if (m_bdivreform)
		{
			if (normE1 > normEm) bconv = false;
		}

		// print convergence summary
		feLog(" Nonlinear solution status: time= %lg\n", tp.currentTime);
		feLog("\tright hand side evaluations   = %d\n", m_nrhs);
		feLog("\tstiffness matrix reformations = %d\n", m_nref);
		feLog("\tconvergence norms :     INITIAL         CURRENT         REQUIRED\n");
		feLog("\t   residual         %15le %15le %15le \n", normRi, normR1, m_Rtol*normRi);
		feLog("\t   energy           %15le %15le %15le \n", normEi, normE1, m_Etol*normEi);
		feLog("\t   displacement     %15le %15le %15le \n", normUi, normu ,(m_Dtol*m_Dtol)*normU );

		// see if we may have a small residual
		if ((bconv == false) && (normR1 < m_Rmin))
		{
			// check for almost zero-residual on the first iteration
			// this might be an indication that there is no force on the system
			feLogWarning("No force acting on the system.");
			bconv = true;
		}

		// see if we have exceeded the max residual
		if ((bconv == false) && (m_Rmax > 0) && (normR1 >= m_Rmax))
		{
			// doesn't look like we're getting anywhere, so let's retry the time step
			throw MaxResidualError();
		}

		// check if we have converged. 
		// If not, calculate the BFGS update vectors
		if (bconv == false)
		{
			if ((normE1 > normEm) && m_bdivreform)
			{
				// check for diverging
				feLogWarning("Problem is diverging. Stiffness matrix will now be reformed");
				normEm = normE1;
				normEi = normE1;
				normRi = normR1;
				m_bforceReform = true;
			}

			// zero displacement increments
			// we must set this to zero before the reformation
			// because we assume that the prescribed displacements are stored 
			// in the m_ui vector.
			zero(m_ui);

			// reform stiffness matrices if necessary
			if (m_bforceReform || m_bdoreforms)
			{
				m_bforceReform = false;

				// reform the matrix
				if (ReformStiffness() == false) return false;
			}

			// copy last calculated residual
			m_R0 = m_R1;
		}
		else if (m_baugment)
		{
			// do the augmentations
			bconv = DoAugmentations();
		}
	
		// increase iteration number
		m_niter++;

		// do minor iterations callbacks
		fem.DoCallback(CB_MINOR_ITERS);
	}
	while (bconv == false);

	// if converged we update the total displacements
	if (bconv)
	{
        UpdateIncrementsEAS(m_Ui, false);
        UpdateIncrements(m_Ut, m_Ui, true);
		zero(m_Ui);
	}

	if (bconv)
	{
		// print a convergence summary to the felog file
		feLog("\nconvergence summary\n");
		feLog("    number of iterations   : %d\n", m_niter);
		feLog("    number of reformations : %d\n", m_nref);
	}

	// if we don't want to hold on to the stiffness matrix, let's clean it up
//	if (m_persistMatrix == false)
	{
		// clean up the solver
		m_plinsolve->Destroy();

		// clean up the stiffness matrix
		m_pK->Clear();

		// make sure we recreate it in the next time step
		m_breshape = true;
	}


	return bconv;
}

//-----------------------------------------------------------------------------
bool FEJacobiSolidSolver::DoAugmentations()
{
	FEModel& fem = *GetFEModel();
	FEAnalysis* pstep = fem.GetCurrentStep();

	// we have converged, so let's see if the augmentations have converged as well
	feLog("\n........................ augmentation # %d\n", m_naug + 1);

	// do callback
	fem.DoCallback(CB_AUGMENT);

	// do the augmentations
	bool bconv = Augment();

	// update counter
	++m_naug;

	// we reset the reformations counter
	m_nref = 0;

	// If we havn't converged we prepare for the next iteration
	if (!bconv)
	{
		// Since the Lagrange multipliers have changed, we can't just copy
		// the last residual but have to recalculate the residual
		// we also recalculate the stresses in case we are doing augmentations
		// for incompressible materials
		fem.Update();
		{
			// TODO: Why is this not calling strategy->Residual? 
			TRACK_TIME(TimerID::Timer_Residual);
			Residual(m_R0);
		}

		// reform the matrix if we are using full-Newton or
		// force reform after augmentations
		if (m_breformAugment)
		{
			// TODO: Note sure how to handle a false return from ReformStiffness. 
			//       I think this is pretty rare so I'm ignoring it for now.
//			if (ReformStiffness() == false) break;
			ReformStiffness();
		}
	}

	return bconv;
}

//-----------------------------------------------------------------------------
//!  Creates the global stiffness matrix
//! \todo Can we move this to the FEGlobalMatrix::Create function?
bool FEJacobiSolidSolver::CreateStiffness(bool breset)
{
	{
		TRACK_TIME(TimerID::Timer_Reform);
		// clean up the solver
		m_plinsolve->Destroy();

		// clean up the stiffness matrix
		m_pK->Clear();

		// create the stiffness matrix
		feLog("===== reforming stiffness matrix:\n");
		if (m_pK->Create(GetFEModel(), m_neq, breset) == false)
		{
			feLogError("An error occured while building the stiffness matrix\n\n");
			return false;
		}
		else
		{
			// output some information about the direct linear solver
			int neq = m_pK->Rows();
			int nnz = m_pK->NonZeroes();
			feLog("\tNr of equations ........................... : %d\n", neq);
			feLog("\tNr of nonzeroes in stiffness matrix ....... : %d\n", nnz);
		}
	}

	// done!
	return true;
}


//-----------------------------------------------------------------------------
//! Reforms a stiffness matrix and factorizes it
bool FEJacobiSolidSolver::ReformStiffness()
{
	feLog("Reforming stiffness matrix: reformation #%d\n\n", m_nref + 1);

	// first, let's make sure we have not reached the max nr of reformations allowed
	if (m_nref >= m_maxref) throw MaxStiffnessReformations();

	FEModel& fem = *GetFEModel();

	// recalculate the shape of the stiffness matrix if necessary
	if (m_breshape)
	{
		// reshape the stiffness matrix
		if (!CreateStiffness(m_niter == 0)) return false;

		// reset reshape flag, except for contact
		m_breshape = (((fem.SurfacePairConstraints() > 0) || (fem.NonlinearConstraints() > 0)) ? true : false);
	}

	// calculate the global stiffness matrix
	bool bret = false;
	{
		TRACK_TIME(TimerID::Timer_Stiffness);

		// zero the stiffness matrix
		m_pK->Zero();

		// Zero the rhs adjustment vector
		zero(m_Fd);

		// calculate the global stiffness matrix
		bret = StiffnessMatrix();

		// check for zero diagonals
		if (m_bzero_diagonal)
		{
			// get the stiffness matrix
			SparseMatrix& K = *m_pK;
			vector<int> zd;
			int neq = K.Rows();
			for (int i = 0; i < neq; ++i)
			{
				double di = fabs(K.diag(i));
				if (di <= m_zero_tol)
				{
					zd.push_back(i);
				}
			}

			if (zd.empty() == false) throw ZeroDiagonal(-1, -1);
		}
	}

	if (bret)
	{
		// increase total nr of reformations
		m_nref++;
		m_ntotref++;
	}

	return bret;
}

//-----------------------------------------------------------------------------
//! Calculates global stiffness matrix.
bool FEJacobiSolidSolver::StiffnessMatrix()
{
	FEModel& fem = *GetFEModel();

	const FETimeInfo& tp = fem.GetTime();

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	// setup the linear system
	FESolidLinearSystem LS(this, &m_rigidSolver, *m_pK, m_Fd, m_ui, (m_msymm == REAL_SYMMETRIC), m_alpha, m_nreq);

	// calculate the stiffness matrix for each domain
	for (int i=0; i<mesh.Domains(); ++i) 
	{
		if (mesh.Domain(i).IsActive()) 
		{
			FEElasticDomain& dom = dynamic_cast<FEElasticDomain&>(mesh.Domain(i));
			dom.StiffnessMatrix(LS);
		}
	}

	// calculate the body force stiffness matrix for each non-rigid domain
	for (int j = 0; j<fem.ModelLoads(); ++j)
	{
		FEModelLoad* pml = fem.ModelLoad(j);
		if (pml->IsActive()) pml->StiffnessMatrix(LS);
	}
    
    // TODO: add body force stiffness for rigid bodies

	// Add mass matrix for dynamic problems
/*	FEAnalysis* pstep = fem.GetCurrentStep();
	if (pstep->m_nanalysis == FESolidAnalysis::DYNAMIC)
	{
		// scale factor
		double dt = tp.timeIncrement;
		double a = tp.alpham / (m_beta*dt*dt);

		// loop over all elastic domains
		for (int i = 0; i<mesh.Domains(); ++i)
		{
			FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
			if (edom) edom->MassMatrix(LS, a);
		}

		m_rigidSolver.RigidMassMatrix(LS, tp);
	}
*/
	// calculate contact stiffness
	ContactStiffness(LS);

	// calculate nonlinear constraint stiffness
	// note that this is the contribution of the 
	// constrainst enforced with augmented lagrangian
	NonLinearConstraintStiffness(LS, tp);

	// add contributions from rigid bodies
	m_rigidSolver.StiffnessMatrix(*m_pK, tp);

	return true;
}

//-----------------------------------------------------------------------------
//! Calculate the stiffness contribution due to nonlinear constraints
void FEJacobiSolidSolver::NonLinearConstraintStiffness(FELinearSystem& LS, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->StiffnessMatrix(LS, tp);
	}
}

//-----------------------------------------------------------------------------
//! This function calculates the contact stiffness matrix

void FEJacobiSolidSolver::ContactStiffness(FELinearSystem& LS)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->StiffnessMatrix(LS, tp);
	}
}

//-----------------------------------------------------------------------------
//! Calculates the contact forces
void FEJacobiSolidSolver::ContactForces(FEGlobalVector& R)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
	{
		FEContactInterface* pci = dynamic_cast<FEContactInterface*>(fem.SurfacePairConstraint(i));
		if (pci->IsActive()) pci->LoadVector(R, tp);
	}
}

//-----------------------------------------------------------------------------
//! calculates the residual vector
//! Note that the concentrated nodal forces are not calculated here.
//! This is because they do not depend on the geometry 
//! so we only calculate them once (in Quasin) and then add them here.

bool FEJacobiSolidSolver::Residual(vector<double>& R)
{
	// get the time information
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();

	// zero nodal reaction forces
	zero(m_Fr);

	// setup the global vector
	zero(R);
	FEResidualVector RHS(fem, R, m_Fr);

	// zero rigid body reaction forces
	m_rigidSolver.Residual();

	// calculate the internal (stress) forces
	InternalForces(RHS);

	// calculate nodal reaction forces
	for (int i = 0; i < m_neq; ++i) m_Fr[i] -= R[i];

	// calculate external forces
	ExternalForces(RHS);

	// increase RHS counter
	m_nrhs++;

	return true;
}

//-----------------------------------------------------------------------------
//! Internal forces
void FEJacobiSolidSolver::InternalForces(FEGlobalVector& R)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	for (int i = 0; i<mesh.Domains(); ++i)
	{
		FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(i));
		if (edom) edom->InternalForces(R);
	}
}

//-----------------------------------------------------------------------------
//! external forces
void FEJacobiSolidSolver::ExternalForces(FEGlobalVector& RHS)
{
	FEModel& fem = *GetFEModel();
	const FETimeInfo& tp = fem.GetTime();
	FEMesh& mesh = fem.GetMesh();

	// apply loads
	for (int j = 0; j<fem.ModelLoads(); ++j)
	{
		FEModelLoad* pml = fem.ModelLoad(j);
		if (pml->IsActive()) pml->LoadVector(RHS);
	}

	// calculate body forces for rigid bodies
	for (int j = 0; j<fem.ModelLoads(); ++j)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.ModelLoad(j));
		if (pbf && pbf->IsActive())
			m_rigidSolver.BodyForces(RHS, tp, *pbf);
	}

	// calculate inertial forces for dynamic problems
	if (fem.GetCurrentStep()->m_nanalysis == FESolidAnalysis::DYNAMIC)
	{
		// allocate F
		vector<double> F;

		// calculate the inertial forces for all elastic domains
		for (int nd = 0; nd < mesh.Domains(); ++nd)
		{
			FEElasticDomain* edom = dynamic_cast<FEElasticDomain*>(&mesh.Domain(nd));
			if (edom) edom->InertialForces(RHS, F);
		}

		// update rigid bodies
		m_rigidSolver.InertialForces(RHS, tp);
	}

	// calculate forces due to surface loads
/*	int nsl = fem.SurfaceLoads();
	for (int i = 0; i<nsl; ++i)
	{
		FESurfaceLoad* psl = fem.SurfaceLoad(i);
		if (psl->IsActive()) psl->LoadVector(RHS, tp);
	}
*/
	// calculate contact forces
	ContactForces(RHS);

	// calculate nonlinear constraint forces
	// note that these are the linear constraints
	// enforced using the augmented lagrangian
	NonLinearConstraintForces(RHS, tp);

	// forces due to point constraints
	//	for (i=0; i<(int) fem.m_PC.size(); ++i) fem.m_PC[i]->Residual(this, R);

	// set the nodal reaction forces
	// TODO: Is this a good place to do this?
	for (int i = 0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		node.set_load(m_dofU[0], 0);
		node.set_load(m_dofU[1], 0);
		node.set_load(m_dofU[2], 0);

		int n;
		if ((n = node.m_ID[m_dofU[0]]) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[0]] - 2) >= 0) node.set_load(m_dofU[0], -m_Fr[n]);

		if ((n = node.m_ID[m_dofU[1]]) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[1]] - 2) >= 0) node.set_load(m_dofU[1], -m_Fr[n]);

		if ((n = node.m_ID[m_dofU[2]]) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);
		if ((n = -node.m_ID[m_dofU[2]] - 2) >= 0) node.set_load(m_dofU[2], -m_Fr[n]);
	}
}

//-----------------------------------------------------------------------------
//! calculate the nonlinear constraint forces 
void FEJacobiSolidSolver::NonLinearConstraintForces(FEGlobalVector& R, const FETimeInfo& tp)
{
	FEModel& fem = *GetFEModel();
	int N = fem.NonlinearConstraints();
	for (int i=0; i<N; ++i) 
	{
		FENLConstraint* plc = fem.NonlinearConstraint(i);
		if (plc->IsActive()) plc->LoadVector(R, tp);
	}
}
