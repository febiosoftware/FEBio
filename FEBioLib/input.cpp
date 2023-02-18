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
#include "FEBioModel.h"
#include <FECore/FEAnalysis.h>
#include <FECore/tens3d.h>
#include <FEBioMech/FERigidJoint.h>
#include <FEBioMech/FERigidSphericalJoint.h>
#include <FEBioPlot/FEBioPlotFile.h>
#include <FECore/FEMaterial.h>
#include <FEBioMech/FERigidBody.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FEInitialCondition.h>
#include <FECore/FENodalLoad.h>
#include <FECore/FESurfaceLoad.h>
#include <FECore/FEBodyLoad.h>
#include <FECore/FEModelParam.h>
#include <FECore/FELoadCurve.h>
#include <FECore/FEMeshAdaptor.h>
#include <FECore/FESurfacePairConstraint.h>
#include <FECore/FETimeStepController.h>
#include <FECore/log.h>
#include <string.h>

//-----------------------------------------------------------------------------
// helper function to print a parameter to the logfile
void FEBioModel::print_parameter(FEParam& p, int level)
{
	char sz[512] = {0};
	int l = (int)strlen(p.name()) + 2*level;
	sprintf(sz, "\t%*s %.*s", l, p.name(), 50-l, "..................................................");
	if (p.dim() == 1)
	{
		switch (p.type())
		{
		case FE_PARAM_DOUBLE : feLog("%s : %lg\n", sz, p.value<double>()); break;
		case FE_PARAM_INT    : 
		{
			if (p.enums())
			{
				int n = p.value<int>();
				const char* szkey = p.enumKey();
				if (szkey)
				{
					feLog("%s : %s (%d)\n", sz, szkey, n);
				}
				else feLog("%s : %d\n", sz, n);
			}
			else feLog("%s : %d\n", sz, p.value<int   >());
		}
		break;
		case FE_PARAM_BOOL   : feLog("%s : %s\n" , sz, (p.value<bool>() ? "yes (1)" : "no (0)")); break;
		case FE_PARAM_STRING : feLog("%s : %s\n" , sz, p.cvalue()); break;
		case FE_PARAM_STD_STRING: feLog("%s : %s\n", sz, p.value<string>().c_str()); break;
		case FE_PARAM_VEC3D  :
			{
				vec3d v = p.value<vec3d>();
				feLog("%s : %lg,%lg,%lg\n", sz, v.x, v.y, v.z);
			}
			break;
		case FE_PARAM_MAT3DS :
			{
				mat3ds m = p.value<mat3ds>();
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg\n", sz, m.xx(), m.yy(), m.zz(), m.xy(), m.yz(), m.xz());
			}
			break;
		case FE_PARAM_MAT3D:
			{
				mat3d m = p.value<mat3d>();
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz, m(0,0), m(0,1), m(0,2), m(1,0), m(1,1), m(1,2), m(2,0), m(2,1), m(2,2));
			}
			break;
		case FE_PARAM_TENS3DRS:
			{
				tens3drs m = p.value<tens3drs>();
				double* d = m.d;
				feLog("%s : %lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg\n", sz,d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10],d[11],d[12],d[13],d[14],d[15],d[16],d[17]);
			}
			break;
		case FE_PARAM_DOUBLE_MAPPED:
			{
				FEParamDouble& v = p.value<FEParamDouble>();
				if (v.isConst())
				{
					double a = v.constValue();
					feLog("%s : %lg\n", sz, a);
				}
				else feLog("%s : (not constant)\n", sz);
			}
			break;
		case FE_PARAM_VEC3D_MAPPED:
		{
			FEParamVec3& v = p.value<FEParamVec3>();
			if (v.isConst())
			{
				vec3d a = v.constValue();
				feLog("%s : %lg, %lg, %lg\n", sz, a.x, a.y, a.z);
			}
			else feLog("%s : (not constant)\n", sz);
		}
		break;
		case FE_PARAM_STD_VECTOR_INT:
		{
			std::vector<int>& v = p.value<std::vector<int> >();
			int nsize = (int)v.size();
			int n = (nsize <= 5 ? nsize : 5);
			feLog("%s : ", sz);
			for (int i = 0; i < n; ++i)
			{
				feLog("%d", v[i]);
				if (i != n - 1) feLog(", ");
			}
			if (nsize > n) feLog(", ...\n");
			else feLog("\n");
		}
		break;
		default:
			feLog("%s : (can't display value)\n", sz);
			break;
		}
	}
	else
	{
		switch (p.type())
		{
		case FE_PARAM_INT:
		case FE_PARAM_DOUBLE:
			{
				int n = p.dim();
				feLog("%s : ", sz);
				for (int k=0; k<n; ++k)
				{
					switch (p.type())
					{
					case FE_PARAM_INT   : feLog("%d" , p.pvalue<int   >()[k]); break;
					case FE_PARAM_DOUBLE: feLog("%lg", p.pvalue<double>()[k]); break;
                    default: break;
					}
					if (k!=n-1) feLog(","); else feLog("\n");
				}
			}
			break;
		case FE_PARAM_DOUBLE_MAPPED:
		{
			int n = p.dim();
			feLog("%s : ", sz);
			for (int k = 0; k < n; ++k)
			{
				FEParamDouble& v = p.value<FEParamDouble>(k);
				if (v.isConst())
					feLog("%lg", v.constValue());
				else
					feLog("???");
				if (k != n - 1) feLog(","); else feLog("\n");
			}
		}
		break;
		default:
			break;
		}
	}
}

//-----------------------------------------------------------------------------
// print the parameter list to the log file
void FEBioModel::print_parameter_list(FEParameterList& pl, int level)
{
	int n = pl.Parameters();
	if (n > 0)
	{
		FEParamIterator it = pl.first();
		for (int j = 0; j < n; ++j, ++it)
		{
			if (it->IsHidden() == false)
			{
				print_parameter(*it, level);
			}
		}
	}
}

//------------------------------------------------------------------------------
void FEBioModel::print_parameter_list(FECoreBase* pc, int level)
{
	print_parameter_list(pc->GetParameterList(), level);
	for (int i=0; i<pc->PropertyClasses(); ++i)
	{
		FEProperty* prop = pc->PropertyClass(i);
		int n = prop->size();
		for (int j=0; j<n; ++j)
		{
			feLog("\t%s: ", prop->GetName());
			FECoreBase* pcj = prop->get(j);
			if (pcj)
			{
				feLog("(type: %s)\n", pcj->GetTypeStr());
				print_parameter_list(pcj, level + 1);
			}
			else feLog("(unspecified)\n");
		}
	}
}

//------------------------------------------------------------------------------
//! This function outputs the input data to the felog file.
void FEBioModel::echo_input()
{
	// get the FE mesh
	FEBioModel& fem = *this;
	FEMesh& mesh = GetMesh();

	// print title
	feLog("%s\n\n", fem.GetTitle().c_str());

	// print file info
	feLog(" FILES USED\n");
	feLog("===========================================================================\n");
	feLog("\tInput file : %s\n", fem.GetInputFileName().c_str());
	feLog("\tPlot file  : %s\n", fem.GetPlotFileName().c_str());
	feLog("\tLog file   : %s\n", fem.GetLogfileName().c_str());
	feLog("\n\n");

	// print mesh info
	feLog(" MESH INFO\n");
	feLog("===========================================================================\n");
	feLog("\tNumber of materials ............................ : %d\n", fem.Materials());
	feLog("\tNumber of domains .............................. : %d\n", mesh.Domains());
	feLog("\tNumber of nodes ................................ : %d\n", mesh.Nodes());
	int nsolid = mesh.Elements(FE_DOMAIN_SOLID    ); if (nsolid > 0) feLog("\tNumber of solid elements ....................... : %d\n", nsolid);
	int nshell = mesh.Elements(FE_DOMAIN_SHELL    ); if (nshell > 0) feLog("\tNumber of shell elements ....................... : %d\n", nshell);
	int nbeam  = mesh.Elements(FE_DOMAIN_BEAM     ); if (nbeam  > 0) feLog("\tNumber of beam elements ........................ : %d\n", nbeam );
	int nelm2d = mesh.Elements(FE_DOMAIN_2D       ); if (nelm2d > 0) feLog("\tNumber of 2D elements .......................... : %d\n", nelm2d);
	feLog("\n\n");

	feLog(" MODULE\n");
	feLog("===========================================================================\n");
	string modName = fem.GetModuleName();
	const char* szmod = modName.c_str();
	if (szmod == 0) { szmod = "unknown"; assert(false); }
	feLog("\tModule type ....................................... : %s\n", szmod);
	feLog("\n\n");

	// print control info
	if (fem.Steps() == 1)
	{
		// get the analysis step
		FEAnalysis& step = *fem.GetStep(0);

		feLog(" CONTROL DATA\n");
		feLog("===========================================================================\n");
		print_parameter_list(step.GetParameterList());

		feLog("\tAuto time stepper activated ....................... : %s\n", (step.m_timeController ? "yes" : "no"));
		if (step.m_timeController)
		{
			FETimeStepController& tc = *step.m_timeController;
			print_parameter_list(tc.GetParameterList(), 1);
		}

		// output solver data
		feLog(" SOLVER PARAMETERS\n");
		feLog("===========================================================================\n");

		FESolver* psolver = step.GetFESolver();
		if (psolver)
		{
			print_parameter_list(psolver->GetParameterList());
			feLog("\n\n");
		}
	}
	else
	{
		feLog(" STEP DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i < fem.Steps(); ++i)
		{
			if (i > 0) feLog("---------------------------------------------------------------------------\n");
			feLog("step %3d - ", i + 1);

			// get the step
			FEAnalysis* pstep = fem.GetStep(i);

			// get the name and type string
			const char* szname = pstep->GetName().c_str();
			const char* sztype = pstep->GetTypeStr();
			if (szname[0] == 0) szname = 0;

			// print type and name
			feLog("%s", (szname ? szname : "(noname)"));
			feLog(" (type: %s)", sztype);
			feLog("\n");

			// print the parameter list
			print_parameter_list(pstep);
		}
	}
	feLog("\n\n");

	// material data
	feLog("\n\n");
	feLog(" MATERIAL DATA\n");
	feLog("===========================================================================\n");
	for (int i=0; i<fem.Materials(); ++i)
	{
		if (i>0) feLog("---------------------------------------------------------------------------\n");
		feLog("%3d - ", i+1);

		// get the material
		FEMaterial* pmat = fem.GetMaterial(i);

		// get the material name and type string
		const char* szname = pmat->GetName().c_str();
		const char* sztype = pmat->GetTypeStr();
		if (szname[0] == 0) szname = 0;

		// print type and name
		feLog("%s", (szname?szname:"unknown"));
		feLog(" (type: %s)", sztype);
		feLog("\n");

		// print the parameter list
		print_parameter_list(pmat);
	}
	feLog("\n\n");

	if (fem.RigidBodies())
	{
		feLog(" RIGID BODY DATA\n");
		feLog("===========================================================================\n");
		for (int i=0; i<fem.RigidBodies(); ++i)
		{
			FERigidBody& rb = *fem.GetRigidBody(i);
			if (i>0) feLog("---------------------------------------------------------------------------\n");
			feLog("Rigid Body %d:\n", rb.m_nID+1);
			feLog("\tmaterial id    : %d\n", rb.m_mat+1);
			feLog("\tname           : %s\n", rb.GetName().c_str());
			feLog("\tcenter of mass : %lg, %lg, %lg\n", rb.m_r0.x, rb.m_r0.y, rb.m_r0.z);
            feLog("\tmass           : %lg\n", rb.m_mass);
            feLog("\tIxx Ixy Ixz    : %lg, %lg, %lg\n", rb.m_moi.xx(), rb.m_moi.xy(), rb.m_moi.xz());
            feLog("\tIxy Iyy Iyz    : %lg, %lg, %lg\n", rb.m_moi.xy(), rb.m_moi.yy(), rb.m_moi.yz());
            feLog("\tIxz Iyz Izz    : %lg, %lg, %lg\n", rb.m_moi.xz(), rb.m_moi.yz(), rb.m_moi.zz());
		}
		feLog("\n\n");
	}

	if (fem.InitialConditions())
	{
		feLog(" INITIAL CONDITION DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i < fem.InitialConditions(); ++i)
		{
			if (i > 0) feLog("---------------------------------------------------------------------------\n");
			feLog("%3d - ", i + 1);

			// get the initial condition
			FEInitialCondition* pic = fem.InitialCondition(i);

			// get the type string
			const char* sztype = pic->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pic->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.BoundaryConditions())
	{
		feLog(" BOUNDARY CONDITION DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i < fem.BoundaryConditions(); ++i)
		{
			if (i > 0) feLog("---------------------------------------------------------------------------\n");
			feLog("%3d - ", i + 1);

			// get the bc
			FEBoundaryCondition* pbc = fem.BoundaryCondition(i);

			// get the type string
			const char* sztype = pbc->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pbc->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.ModelLoads())
	{
		feLog(" MODEL LOAD DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i < fem.ModelLoads(); ++i)
		{
			if (i > 0) feLog("---------------------------------------------------------------------------\n");
			feLog("%3d - ", i + 1);

			// get the model load
			FEModelLoad* pml = fem.ModelLoad(i);

			// get the type string
			const char* sztype = pml->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog(" Type: %s\n", sztype);

			// print the parameter list
			FEParameterList& pl = pml->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.SurfacePairConstraints() > 0)
	{
		feLog(" CONTACT INTERFACE DATA\n");
		feLog("===========================================================================\n");
		for (int i = 0; i<fem.SurfacePairConstraints(); ++i)
		{
			if (i>0) feLog("---------------------------------------------------------------------------\n");

			FESurfacePairConstraint* pi = fem.SurfacePairConstraint(i);
			const char* sztype = pi->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog("contact interface %d - Type: %s\n", i+1, sztype);
			FEParameterList& pl = pi->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.NonlinearConstraints() != 0)
	{
		feLog(" NONLINEAR CONSTRAINT DATA\n");
		feLog("===========================================================================\n");
		int NC = fem.NonlinearConstraints();
		for (int i = 0; i<NC; ++i)
		{
			FENLConstraint* plc = fem.NonlinearConstraint(i); assert(plc);
			const char* sztype = plc->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog("\nnonlinear constraint %d - Type: %s\n", i + 1, sztype);
			FEParameterList& pl = plc->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	if (fem.MeshAdaptors())
	{
		feLog(" MESH ADAPTOR DATA\n");
		feLog("===========================================================================\n");
		int NMA = fem.MeshAdaptors();
		for (int i = 0; i < NMA; ++i)
		{
			if (i > 0) feLog("---------------------------------------------------------------------------\n");
			FEMeshAdaptor* pma = fem.MeshAdaptor(i); assert(pma);
			const char* sztype = pma->GetTypeStr();
			if (sztype == 0) sztype = "unknown";
			feLog("\nmesh adaptor %d - Type: %s\n", i + 1, sztype);
			FEParameterList& pl = pma->GetParameterList();
			print_parameter_list(pl);
		}
		feLog("\n\n");
	}

	feLog(" LOAD CONTROLLER DATA\n");
	feLog("===========================================================================\n");
	for (int i = 0; i < fem.LoadControllers(); ++i)
	{
		if (i > 0) feLog("---------------------------------------------------------------------------\n");
		feLog("load controller %3d - ", i + 1);

		// get the load controller
		FELoadController* plc = fem.GetLoadController(i);

		// get the name and type string
		const char* szname = plc->GetName().c_str();
		const char* sztype = plc->GetTypeStr();
		if (szname[0] == 0) szname = 0;

		// print type and name
		feLog("%s", (szname ? szname : "(noname)"));
		feLog(" (type: %s)", sztype);
		feLog("\n");

		// print the parameter list
		print_parameter_list(plc);
	}
	feLog("\n\n");

	// print output data
	feLog(" OUTPUT DATA\n");
	feLog("===========================================================================\n");

	PlotFile* pplt = fem.GetPlotFile();
	if (dynamic_cast<FEBioPlotFile*>(pplt))
	{
		FEBioPlotFile* pf = dynamic_cast<FEBioPlotFile*>(pplt);
		feLog("\tplotfile format ........................... : FEBIO\n");

		const FEBioPlotFile::Dictionary& dic = pf->GetDictionary();

		for (int i = 0; i < 5; ++i)
		{
			const list<FEBioPlotFile::DICTIONARY_ITEM>* pl = 0;
			const char* szn = 0;
			switch (i)
			{
			case 0: pl = &dic.GlobalVariableList(); szn = "Global Variables"; break;
			case 1: pl = &dic.MaterialVariableList(); szn = "Material Variables"; break;
			case 2: pl = &dic.NodalVariableList(); szn = "Nodal Variables"; break;
			case 3: pl = &dic.DomainVariableList(); szn = "Domain Variables"; break;
			case 4: pl = &dic.SurfaceVariableList(); szn = "Surface Variables"; break;
			}

			if (!pl->empty())
			{
				feLog("\t\t%s:\n", szn);
				list<FEBioPlotFile::DICTIONARY_ITEM>::const_iterator it;
				for (it = pl->begin(); it != pl->end(); ++it)
				{
					const char* szt = 0;
					switch (it->m_ntype)
					{
					case PLT_FLOAT: szt = "float"; break;
					case PLT_VEC3F: szt = "vec3f"; break;
					case PLT_MAT3FS: szt = "mat3fs"; break;
					case PLT_MAT3FD: szt = "mat3fd"; break;
					case PLT_TENS4FS: szt = "tens4fs"; break;
					case PLT_MAT3F: szt = "mat3f"; break;
					}

					const char* szf = 0;
					switch (it->m_nfmt)
					{
					case FMT_NODE: szf = "NODE"; break;
					case FMT_ITEM: szf = "ITEM"; break;
					case FMT_MULT: szf = "COMP"; break;
					case FMT_REGION: szf = "REGION"; break;
					}

					feLog("\t\t\t%-20s (type = %5s, format = %4s)\n", it->m_szname, szt, szf);
				}
			}
		}
	}

	FECoreKernel& fecore = FECoreKernel::GetInstance();
	feLog(" LINEAR SOLVER DATA\n");
	feLog("===========================================================================\n");
	feLog("\tDefault linear solver ............................. : %s\n", fecore.GetLinearSolverType());
	FECoreFactory* fac = fecore.FindFactoryClass(FELINEARSOLVER_ID, fecore.GetLinearSolverType());
	if (fac) print_parameter_list(fac->GetParameterList());
	feLog("\n");
}
