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
#include "FEBioRVE.h"
#include "FEMicroMaterial.h"
#include "FEMicroMaterial2O.h"
#include "FEMindlinElastic2O.h"
#include "FEPeriodicBoundary2O.h"
#include "FE2OMicroConstraint.h"
#include "FEElasticMultiscaleDomain1O.h"
#include "FEElasticMultiscaleDomain2O.h"
#include "FEMultiscaleDomainFactory.h"
#include "FEBioRVEPlot.h"
#include "FERVEProbe.h"

//-----------------------------------------------------------------------------
//! Register all the classes of the FEBioMech module with the FEBio framework.
void FEBioRVE::InitModule()
{
	FECoreKernel& febio = FECoreKernel::GetInstance();
	febio.RegisterDomain(new FEMultiScaleDomainFactory, true);

	// this module extends the solid module
	febio.SetActiveModule("solid");

	REGISTER_FECORE_CLASS(FEMicroMaterial, "micro-material");
	REGISTER_FECORE_CLASS(FEMicroMaterial2O, "micro-material2O");
	REGISTER_FECORE_CLASS(FEMindlinElastic2O, "mindlin elastic");

	REGISTER_FECORE_CLASS(FEMicroProbe, "probe");

	REGISTER_FECORE_CLASS(FEElasticMultiscaleDomain1O, "elastic-mm-solid");
	REGISTER_FECORE_CLASS(FEElasticMultiscaleDomain2O, "elastic-mm-solid2O");
	REGISTER_FECORE_CLASS(FEElasticSolidDomain2O, "elastic-solid2O");

	REGISTER_FECORE_CLASS(FE2OMicroConstraint, "2O microfluc");

	REGISTER_FECORE_CLASS(FEPeriodicBoundary1O, "periodic boundary1O");
	REGISTER_FECORE_CLASS(FEPeriodicBoundary2O, "periodic boundary2O");

	REGISTER_FECORE_CLASS(FEPlotElementGnorm, "G norm");
	REGISTER_FECORE_CLASS(FEPlotElementPK1norm, "PK1 norm");
	REGISTER_FECORE_CLASS(FEPlotElementQK1norm, "QK1 norm");
	REGISTER_FECORE_CLASS(FEPlotElementMicroEnergy, "micro energy");
}
