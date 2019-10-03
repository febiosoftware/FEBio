/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2019 University of Utah, The Trustees of Columbia University in
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
#include "FEFluidSolutesSolver2.h"
#include "FEFluid.h"
#include "FESolutesMaterial.h"
#include "FEFluidSolutesDomain2.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include <FECore/log.h>


BEGIN_FECORE_CLASS(FEFluidSolutesSolver2, FESolver)
	ADD_PARAMETER(m_maxRef, "max_refs");
	ADD_PARAMETER(m_maxUps, "max_ups");
	ADD_PARAMETER(m_divergeReform, "diverge_reform");
	ADD_PARAMETER(m_reformTimeStep, "reform_each_time_step");
	ADD_PARAMETER(m_vtol, "vtol");
	ADD_PARAMETER(m_ftol, "ftol");
	ADD_PARAMETER(m_ctol, "ctol");
	ADD_PARAMETER(m_etol, "etol");
	ADD_PARAMETER(m_rtol, "rtol");
	ADD_PARAMETER(m_lstol, "lstol");
	ADD_PARAMETER(m_minRes, "min_residual");
	ADD_PARAMETER(m_rhoi, "rhoi");
	ADD_PARAMETER(m_qnmethod, "qnmethod");
END_FECORE_CLASS();

//! constructor
FEFluidSolutesSolver2::FEFluidSolutesSolver2(FEModel* fem) : FESolver(fem), m_fldSolver(fem), m_sltSolver(fem)
{
	m_maxRef = 15;
	m_maxUps = 10;
	m_divergeReform = true;
	m_reformTimeStep = true;
	m_etol = 0.01;
	m_rtol = 0.0;
	m_lstol = 0.0;
	m_minRes = 0.0;
	m_rhoi = 0.0;
	m_qnmethod = QN_BROYDEN;

	m_vtol = 0.01;
	m_ftol = 0.01;
	m_ctol = 0.01;
}

//! destructor
FEFluidSolutesSolver2::~FEFluidSolutesSolver2()
{

}

//! Initializes data structures
bool FEFluidSolutesSolver2::Init()
{
	if (FESolver::Init() == false) return false;

	// assign fluid solver parameters
	m_fldSolver.m_maxref = m_maxRef;
	m_fldSolver.m_maxups = m_maxUps;
	m_fldSolver.m_bdivreform = m_divergeReform;
	m_fldSolver.m_breformtimestep = m_reformTimeStep;
	m_fldSolver.m_Vtol = m_vtol;
	m_fldSolver.m_Ftol = m_ftol;
	m_fldSolver.m_Etol = m_etol;
	m_fldSolver.m_Rtol = m_rtol;
	m_fldSolver.m_lineSearch->m_LStol = m_lstol;
	m_fldSolver.m_Rmin = m_minRes;
	m_fldSolver.m_rhoi = m_rhoi;
	m_fldSolver.m_qndefault = m_qnmethod;

	// assign solute solver parameters
	m_sltSolver.m_maxref = m_maxRef;
	m_sltSolver.m_maxups = m_maxUps;
	m_sltSolver.m_bdivreform = m_divergeReform;
	m_sltSolver.m_breformtimestep = m_reformTimeStep;
	m_sltSolver.m_Ctol = m_ctol;
	m_sltSolver.m_Etol = m_etol;
	m_sltSolver.m_Rtol = m_rtol;
	m_sltSolver.m_lineSearch->m_LStol = m_lstol;
	m_sltSolver.m_Rmin = m_minRes;
	m_sltSolver.m_rhoi = m_rhoi;
	m_sltSolver.m_qndefault = m_qnmethod;

	// We need to enforce the block allocation strategy
	m_fldSolver.SetEquationScheme(EQUATION_SCHEME::BLOCK);
	m_sltSolver.SetEquationScheme(EQUATION_SCHEME::BLOCK);

	// initialize fluid solver
	if (m_fldSolver.Init() == false) return false;

	// initialize solute solver
	if (m_sltSolver.Init() == false) return false;

	return true;
}

//! initialize the step
bool FEFluidSolutesSolver2::InitStep(double time)
{
	// initialize fluid solver
	if (m_fldSolver.InitStep(time) == false) return false;

	// initialize solute solver
	if (m_sltSolver.InitStep(time) == false) return false;

	return true;
}

//! Initialize linear equation system
bool FEFluidSolutesSolver2::InitEquations()
{
	// TODO: I think this is going to call FENewtonSolver::InitEquations twice
	// but I don't think this is an issue since the same equations will be allocated
	// We do need to call the InitEquations of each solver to initialize local
	// equation counters.

	// initialize fluid solver
	if (m_fldSolver.InitEquations2() == false) return false;

	// initialize solute solver
	if (m_sltSolver.InitEquations2() == false) return false;
	
	return true;
}

bool FEFluidSolutesSolver2::SolveStep()
{
	// get the mesh
	FEMesh& mesh = GetFEModel()->GetMesh();

	// active the fluid domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEFluidSolutesDomain2& dom = dynamic_cast<FEFluidSolutesDomain2&>(mesh.Domain(i));
		dom.SetActiveDomain(FEFluidSolutesDomain2::FLUID_DOMAIN);
	}

	// First, the fluid equations are solved
	feLog("\n** Solving for velocity field:\n");
	if (m_fldSolver.SolveStep() == false) return false;

	// next, we map the fluid velocity and its divergence to the integration points
	// of the solute model
	MapVelocitySolution();

	// active the solute domains
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEFluidSolutesDomain2& dom = dynamic_cast<FEFluidSolutesDomain2&>(mesh.Domain(i));
		dom.SetActiveDomain(FEFluidSolutesDomain2::SOLUTES_DOMAIN);
	}

	// Now, we solve the solutes problem
	feLog("\n** Solving for solute concentrations:\n");
	if (m_sltSolver.SolveStep() == false) return false;

	// All good!
	return true;
}

//! Serialization
void FEFluidSolutesSolver2::Serialize(DumpStream& ar)
{
	m_fldSolver.Serialize(ar);
	m_sltSolver.Serialize(ar);
}

void FEFluidSolutesSolver2::MapVelocitySolution()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// loop over all domains
	for (int m = 0; m < mesh.Domains(); ++m)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(m));

		const int NE = dom.Elements();
		for (int i = 0; i < NE; ++i)
		{
			FESolidElement& el = dom.Element(i);

			const int nint = el.GaussPoints();
			for (int n = 0; n < nint; ++n)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(n);
				FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
				FESolutesMaterial::Point& sp = *mp.ExtractData<FESolutesMaterial::Point>();

				sp.m_vft = fp.m_vft;
				sp.m_divf = fp.m_Jfdot / fp.m_Jf;
			}			
		}
	}
}
