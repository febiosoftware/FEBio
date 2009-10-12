// FEBioImport.cpp: implementation of the FEFEBioImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FEBioImport.h"
#include "FERigid.h"
#include "FEFacet2FacetSliding.h"
#include "FESlidingInterface2.h"
#include "FECore/ConjGradIterSolver.h"
#include "FESolidSolver.h"
#include "FEHeatSolver.h"
#include <string.h>

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FEFEBioImport::Load
//  Imports an XML input file.
//  The actual file is parsed using the XMLReader class.
//

bool FEFEBioImport::Load(FEM& fem, const char* szfile)
{
	// store a copy of the file name
	fem.SetInputFilename(szfile);

	// Open the XML file
	if (m_xml.Open(szfile) == false) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// keep a pointer to the fem object
	m_pfem = &fem;

	// get a pointer to the first step
	// since we assume that the FEM object will always have
	// at least one step defined
	int nsteps = fem.m_Step.size();
	assert(nsteps > 0);
	m_pStep = &fem.m_Step[nsteps-1];
	m_nsteps = 0; // reset step section counter

	// get the logfile
	Logfile& log = fem.m_log;

	// loop over all child tags
	try
	{
		// Find the root element
		XMLTag tag;
		if (m_xml.FindTag("febio_spec", tag) == false) return false;

		if (strcmp(tag.m_szatv[0], "1.0") == 0)
		{
			// Read version 1.0

			++tag;
			do
			{
				if		(tag == "Module"  ) ParseModuleSection  (tag);
				else if (tag == "Control" ) ParseControlSection (tag);
				else if (tag == "Material") ParseMaterialSection(tag);
				else if (tag == "Geometry") ParseGeometrySection(tag);
				else if (tag == "Boundary") ParseBoundarySection(tag);
				else if (tag == "Initial" ) ParseInitialSection (tag);
				else if (tag == "Globals" ) ParseGlobalsSection (tag);
				else if (tag == "LoadData") ParseLoadSection    (tag);
				else if (tag == "Output"  ) ParseOutputSection  (tag);
				else if (tag == "Step"    ) ParseStepSection    (tag);
				else throw XMLReader::InvalidTag(tag);

				// go to the next tag
				++tag;
			}
			while (!tag.isend());
		}
		else throw InvalidVersion();
	}
	catch (InvalidVersion)
	{
		log.printbox("FATAL ERROR", "Invalid version for FEBio specification.");
		return false;
	}
	catch (InvalidMaterial e)
	{
		log.printbox("FATAL ERROR:", "Element %d has an invalid material type.", e.m_nel);
		return false;
	}
	catch (XMLReader::XMLSyntaxError)
	{
		log.printf("FATAL ERROR: Syntax error (line %d)\n", m_xml.GetCurrentLine());
		return false;
	}
	catch (XMLReader::InvalidTag e)
	{
		log.printf("FATAL ERROR: unrecognized tag \"%s\" (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::InvalidAttributeValue e)
	{
		const char* szt = e.tag.m_sztag;
		const char* sza = e.szatt;
		const char* szv = e.szval;
		int l = e.tag.m_nstart_line;
		log.printf("FATAL ERROR: unrecognized value \"%s\" for attribute \"%s.%s\" (line %d)\n", szv, szt, sza, l);
		return false;
	}
	catch (XMLReader::InvalidValue e)
	{
		log.printf("FATAL ERROR: the value for tag \"%s\" is invalid (line %d)\n", e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::MissingAttribute e)
	{
		log.printf("FATAL ERROR: Missing attribute \"%s\" of tag \"%s\" (line %d)\n", e.szatt, e.tag.m_sztag, e.tag.m_nstart_line);
		return false;
	}
	catch (XMLReader::UnmatchedEndTag e)
	{
		const char* sz = e.tag.m_szroot[e.tag.m_nlevel];
		log.printf("FATAL ERROR: Unmatched end tag for \"%s\" (line %d)\n", sz, e.tag.m_nstart_line);
		return false;
	}
	catch (...)
	{
		log.printf("FATAL ERROR: unrecoverable error (line %d)\n", m_xml.GetCurrentLine());
		return false;
	}

	// close the XML file
	m_xml.Close();

	// we're done!
	return true;
}

//-----------------------------------------------------------------------------
//! This function parses the Module section.
//! The Module defines the type of problem the user wants to solve (solid, heat, ...)
//!

bool FEFEBioImport::ParseModuleSection(XMLTag &tag)
{
	FEM& fem = *m_pfem;

	// get the type attribute
	const char* szt = tag.AttributeValue("type");

	assert(m_pStep && (m_pStep->m_psolver == 0));

	if (strcmp(szt, "solid") == 0) m_pStep->m_psolver = new FESolidSolver(fem);
	else if (strcmp(szt, "heat") == 0)
	{
		m_pStep->m_psolver = new FEHeatSolver(fem);
		m_pStep->m_nModule = FE_HEAT;
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseControlSection
//  This function parses the control section from the xml file
//

bool FEFEBioImport::ParseControlSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;

	// make sure we have a solver defined
	if (m_pStep->m_psolver == 0) m_pStep->m_psolver = new FESolidSolver(fem);

	char sztitle[256];

	++tag;
	do
	{
		if      (tag == "title"             ) { tag.value(sztitle); fem.SetTitle(sztitle); }
		else if (tag == "time_steps"        ) tag.value(m_pStep->m_ntime);
		else if (tag == "step_size"         ) tag.value(m_pStep->m_dt0);
		else if (tag == "dtol"              ) tag.value(m_pStep->m_psolver->m_Dtol);
		else if (tag == "etol"              ) tag.value(m_pStep->m_psolver->m_Etol);
		else if (tag == "rtol"              ) tag.value(m_pStep->m_psolver->m_Rtol);
		else if (tag == "lstol"             ) tag.value(m_pStep->m_psolver->m_LStol);
		else if (tag == "lsmin"             ) tag.value(m_pStep->m_psolver->m_LSmin);
		else if (tag == "lsiter"            ) tag.value(m_pStep->m_psolver->m_LSiter);
		else if (tag == "max_refs"          ) tag.value(m_pStep->m_psolver->m_maxref);
		else if (tag == "max_ups"           ) tag.value(m_pStep->m_psolver->m_maxups);
		else if (tag == "cmax"              ) tag.value(m_pStep->m_psolver->m_cmax);
		else if (tag == "optimize_bw"       ) tag.value(fem.m_bwopt);
		else if (tag == "pressure_stiffness") tag.value(m_pStep->m_istiffpr);
		else if (tag == "hourglass"         ) tag.value(m_pStep->m_hg);
		else if (tag == "symmetric_biphasic") tag.value(fem.m_bsym_poro);
		else if (tag == "plane_strain"      ) tag.value(fem.m_bplane_strain);
		else if (tag == "analysis")
		{
			const char* szt = tag.AttributeValue("type");
			if      (strcmp(szt, "static" ) == 0) m_pStep->m_nanalysis = FE_STATIC;
			else if (strcmp(szt, "dynamic") == 0) m_pStep->m_nanalysis = FE_DYNAMIC;
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
		else if (tag == "restart" )
		{
			const char* szf = tag.AttributeValue("file", true);
			if (szf) fem.SetDumpFilename(szf);
			tag.value(m_pStep->m_bDump);
		}
		else if (tag == "time_stepper")
		{
			m_pStep->m_bautostep = true;
			++tag;
			do
			{
				if      (tag == "max_retries") tag.value(m_pStep->m_maxretries);
				else if (tag == "opt_iter"   ) tag.value(m_pStep->m_iteopt);
				else if (tag == "dtmin"      ) tag.value(m_pStep->m_dtmin);
				else if (tag == "dtmax"      )
				{
					tag.value(m_pStep->m_dtmax);
					const char* sz = tag.AttributeValue("lc", true);
					if (sz) sscanf(sz,"%d", &m_pStep->m_nmplc);
				}
				else if (tag == "aggressiveness") tag.value(m_pStep->m_naggr);
				else throw XMLReader::InvalidTag(tag);

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "plot_level")
		{
			char szval[256];
			tag.value(szval);
			if		(strcmp(szval, "PLOT_DEFAULT"    ) == 0) {} // don't change the plot level
			else if (strcmp(szval, "PLOT_NEVER"      ) == 0) m_pStep->SetPlotLevel(FE_PLOT_NEVER);
			else if (strcmp(szval, "PLOT_MAJOR_ITRS" ) == 0) m_pStep->SetPlotLevel(FE_PLOT_MAJOR_ITRS);
			else if (strcmp(szval, "PLOT_MINOR_ITRS" ) == 0) m_pStep->SetPlotLevel(FE_PLOT_MINOR_ITRS);
			else if (strcmp(szval, "PLOT_MUST_POINTS") == 0) m_pStep->SetPlotLevel(FE_PLOT_MUST_POINTS);
			else throw XMLReader::InvalidValue(tag);
		}
		else if (tag == "print_level")
		{
			char szval[256];
			tag.value(szval);
			if      (strcmp(szval, "PRINT_DEFAULT"       ) == 0) {} // don't change the default print level
			else if (strcmp(szval, "PRINT_NEVER"         ) == 0) m_pStep->SetPrintLevel(FE_PRINT_NEVER);
			else if (strcmp(szval, "PRINT_PROGRESS"      ) == 0) m_pStep->SetPrintLevel(FE_PRINT_PROGRESS);
			else if (strcmp(szval, "PRINT_MAJOR_ITRS"    ) == 0) m_pStep->SetPrintLevel(FE_PRINT_MAJOR_ITRS);
			else if (strcmp(szval, "PRINT_MINOR_ITRS"    ) == 0) m_pStep->SetPrintLevel(FE_PRINT_MINOR_ITRS);
			else if (strcmp(szval, "PRINT_MINOR_ITRS_EXP") == 0) m_pStep->SetPrintLevel(FE_PRINT_MINOR_ITRS_EXP);
			else throw XMLReader::InvalidTag(tag);
		}
		else if (tag == "integration")
		{
			++tag;
			do
			{
				if (tag == "rule")
				{
					const char* sze = tag.AttributeValue("elem");
					const char* szv = tag.szvalue();

					if (strcmp(sze, "hex8") == 0)
					{
						if (strcmp(szv, "GAUSS8") == 0) fem.m_nhex8 = FE_HEX;
						else if (strcmp(szv, "POINT6") == 0) fem.m_nhex8 = FE_RIHEX;
						else if (strcmp(szv, "UDG") == 0) fem.m_nhex8 = FE_UDGHEX;
						else throw XMLReader::InvalidValue(tag);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "elem", sze);
				}
				else throw XMLReader::InvalidValue(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "linear_solver")
		{
			const char* szt = tag.AttributeValue("type");
			if      (strcmp(szt, "skyline"           ) == 0) fem.m_nsolver = SKYLINE_SOLVER;
			else if (strcmp(szt, "psldlt"            ) == 0) fem.m_nsolver = PSLDLT_SOLVER;
			else if (strcmp(szt, "superlu"           ) == 0) fem.m_nsolver = SUPERLU_SOLVER;
			else if (strcmp(szt, "superlu_mt"        ) == 0) fem.m_nsolver = SUPERLU_MT_SOLVER;
			else if (strcmp(szt, "pardiso"           ) == 0) fem.m_nsolver = PARDISO_SOLVER;
			else if (strcmp(szt, "wsmp"              ) == 0) fem.m_nsolver = WSMP_SOLVER;
			else if (strcmp(szt, "lusolver"          ) == 0) fem.m_nsolver = LU_SOLVER;
			else if (strcmp(szt, "conjugate gradient") == 0)
			{
				fem.m_nsolver = CG_ITERATIVE_SOLVER;
				ConjGradIterSolver* ps;
				fem.m_pStep->m_psolver->m_psolver = ps = new ConjGradIterSolver();
				if (!tag.isleaf())
				{
					++tag;
					do
					{
						if      (tag == "tolerance"     ) tag.value(ps->m_tol);
						else if (tag == "max_iterations") tag.value(ps->m_kmax);
						else if (tag == "print_level"   ) tag.value(ps->m_nprint);
						else throw XMLReader::InvalidTag(tag);
						++tag;
					}
					while (!tag.isend());
				}
			}
			else throw XMLReader::InvalidAttributeValue(tag, "type", szt);
		}
		else throw XMLReader::InvalidTag(tag);

		if (fem.m_pStep->m_psolver->m_Rtol == 0) fem.m_pStep->m_psolver->m_Rtol = 1e10;

		++tag;
	}
	while (!tag.isend());

	// add the "zero" loadcurve
	// this is the loadcurve that will be used if a loadcurve is not
	// specified for something that depends on time
	FELoadCurve* plc = new FELoadCurve;

	plc->Create(2);
	plc->LoadPoint(0).time = 0;
	plc->LoadPoint(0).value = 0;
	plc->LoadPoint(1).time = m_pStep->m_ntime*m_pStep->m_dt0;
	plc->LoadPoint(1).value = 1;

	fem.AddLoadCurve(plc);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseMaterialSection
//  This function parses the material section from the xml file
//

bool FEFEBioImport::ParseMaterialSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;

	const char* sztype = 0;
	const char* szname = 0;

	int nmat = 0;

	++tag;
	do
	{
		// get the material type
		sztype = tag.AttributeValue("type");

		// get the material name
		szname = tag.AttributeValue("name", true);

		// create a new material of this type
		FEMaterial* pmat = FEMaterialFactory::CreateMaterial(sztype);
		if (pmat == 0) throw XMLReader::InvalidAttributeValue(tag, "type", sztype);

		// add the material
		fem.AddMaterial(pmat);
		++nmat;

		// set the material's name
		if (szname) pmat->SetName(szname);

		// get the material's parameter list
		FEParameterList* pl = pmat->GetParameterList();
		fem.AddParameterList(pl);

		// loop over all parameters
		++tag;
		do
		{
			// see if we can find this parameter
			FEParam* pp = pl->Find(tag.Name());
			if (pp)
			{
				switch (pp->m_itype)
				{
				case FE_PARAM_DOUBLE : tag.value(pp->value<double>() ); break;
				case FE_PARAM_INT    : tag.value(pp->value<int   >() ); break;
				case FE_PARAM_BOOL   : tag.value(pp->value<bool  >() ); break;
				case FE_PARAM_INTV   : tag.value(pp->pvalue<int   >(), pp->m_ndim); break;
				case FE_PARAM_DOUBLEV: tag.value(pp->pvalue<double>(), pp->m_ndim); break;
				default:
					assert(false);
				}

				int lc = -1;
				tag.AttributeValue("lc", lc, true);
				if (lc != -1) pp->m_nlc = lc;
			}
			else
			{
				// if we get here the parameter was not part of the parameter list
				// however, not all parameters can be read from the parameter lists yet
				// so we have to read in the other parameters the hard way
				// TODO: this option will become obselete
				bool bfound = false;

				// additional elastic material parameters
				if (!bfound && dynamic_cast<FEElasticMaterial*>(pmat))
				{
					FEElasticMaterial* pm = dynamic_cast<FEElasticMaterial*>(pmat);

					// read the material axis
					if (tag == "mat_axis")
					{
						const char* szt = tag.AttributeValue("type");
						if (strcmp(szt, "local") == 0)
						{
							FELocalMap* pmap = new FELocalMap();
							pm->m_pmap = pmap;

							int n[3];
							tag.value(n, 3);
							if ((n[0] == 0) && (n[1] == 0) && (n[2] == 0)) { n[0] = 1; n[1] = 2; n[2] = 4; }

							pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
						}
						else if (strcmp(szt, "vector") == 0)
						{
							FEVectorMap* pmap = new FEVectorMap();
							pm->m_pmap = pmap;

							vec3d a(1,0,0), d(0,1,0);
							++tag;
							do
							{
								if (tag == "a") tag.value(a);
								else if (tag == "d") tag.value(d);
								else throw XMLReader::InvalidTag(tag);

								++tag;
							}
							while (!tag.isend());
							pmap->SetVectors(a, d);
						}
						else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

						bfound = true;
					}
				}

				// additional transversely isotropic material parameters
				if (!bfound && dynamic_cast<FETransverselyIsotropic*>(pmat))
				{
					// read the transversely isotropic materials data
					FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*>(pmat);

					// read material fibers
					if (tag == "fiber")
					{
						const char* szt = tag.AttributeValue("type");
						if (strcmp(szt, "local") == 0)
						{
							FELocalMap* pmap = new FELocalMap();
							pm->m_pmap = pmap;

							int n[3] = {0};
							tag.value(n, 2);

							if ((n[0]==0)&&(n[1]==0)&&(n[2]==0)) { n[0] = 1; n[1] = 2; n[2] = 4; }
							if (n[2] == 0) n[2] = n[1];

							pmap->SetLocalNodes(n[0]-1, n[1]-1, n[2]-1);
						}
						else if (strcmp(szt, "spherical") == 0)
						{
							FESphericalMap* pmap = new FESphericalMap();
							pm->m_pmap = pmap;

							vec3d c;
							tag.value(c);

							pmap->SetSphereCenter(c);
						}
						else if (strcmp(szt, "vector") == 0)
						{
							FEVectorMap* pmap = new FEVectorMap;
							pm->m_pmap = pmap;

							vec3d a, d;
							tag.value(a);
							a.unit();

							d = vec3d(1,0,0);
							if (a*d > .999) d = vec3d(0,1,0);

							pmap->SetVectors(a, d);
						}
						else if (strcmp(szt, "random2d") == 0)
						{
							FERandom2DMap* pmap = new FERandom2DMap();
							pm->m_pmap = pmap;
						}
						else if (strcmp(szt, "user") == 0)
						{
							// fibers are read in in the ElementData section
						}
						else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

						// mark the tag as read
						bfound = true;
					}
					else if (tag == "active_contraction")
					{
						const char* szlc = tag.AttributeValue("lc", true);
						int lc = 0;
						if (szlc) lc = atoi(szlc);
						pm->lcna = lc;
						tag.value(pm->m_ascl);

						if (!tag.isleaf())
						{
							++tag;
							do
							{
								if (tag == "ca0") tag.value(pm->ca0);
								else if (tag == "beta") tag.value(pm->beta);
								else if (tag == "l0") tag.value(pm->l0);
								else if (tag == "refl") tag.value(pm->refl);
								else throw XMLReader::InvalidTag(tag);
								++tag;
							}
							while (!tag.isend());
						}

						// mark tag as read
						bfound = true;
					}
				}

				// read rigid body data
				if (!bfound && dynamic_cast<FERigid*>(pmat))
				{
					FERigid* pm = dynamic_cast<FERigid*>(pmat);

					if (tag == "center_of_mass") { tag.value(pm->m_rc); pm->m_com = 1; bfound = true; }
					else if (strncmp(tag.Name(), "trans_", 6) == 0)
					{
						const char* szt = tag.AttributeValue("type");
						const char* szlc = tag.AttributeValue("lc", true);
						int lc = 0;
						if (szlc) lc = atoi(szlc)+1;

						int bc = -1;
						if      (tag.Name()[6] == 'x') bc = 0;
						else if (tag.Name()[6] == 'y') bc = 1;
						else if (tag.Name()[6] == 'z') bc = 2;
						assert(bc >= 0);

						if      (strcmp(szt, "free"      ) == 0) pm->m_bc[bc] =  0;
						else if (strcmp(szt, "fixed"     ) == 0) pm->m_bc[bc] = -1;
						else if (strcmp(szt, "prescribed") == 0)
						{
							pm->m_bc[bc] = lc;
							FERigidBodyDisplacement DC;
							DC.id = nmat;
							DC.bc = bc;
							DC.lc = lc;
							tag.value(DC.sf);
							fem.m_RDC.add(DC);

							// add this boundary condition to the current step
							if (m_nsteps > 0)
							{
								int n = fem.m_RDC.size()-1;
								FERigidBodyDisplacement* pDC = &fem.m_RDC[n];
								m_pStep->AddBoundaryCondition(pDC);
								pDC->Deactivate();
							}
						}
						else if (strcmp(szt, "force") == 0)
						{
							pm->m_bc[bc] = 0;
							FERigidBodyForce FC;
							FC.id = nmat;
							FC.bc = bc;
							FC.lc = lc-1;
							tag.value(FC.sf);
							fem.m_RFC.add(FC);

							// add this boundary condition to the current step
							if (m_nsteps > 0)
							{
								int n = fem.m_RFC.size()-1;
								FERigidBodyForce* pFC = &fem.m_RFC[n];
								m_pStep->AddBoundaryCondition(pFC);
								pFC->Deactivate();
							}
						}
						bfound = true;
					}
					else if (strncmp(tag.Name(), "rot_", 4) == 0)
					{
						const char* szt = tag.AttributeValue("type");
						const char* szlc = tag.AttributeValue("lc", true);
						int lc = 0;
						if (szlc) lc = atoi(szlc)+1;

						int bc = -1;
						if      (tag.Name()[4] == 'x') bc = 3;
						else if (tag.Name()[4] == 'y') bc = 4;
						else if (tag.Name()[4] == 'z') bc = 5;
						assert(bc >= 0);

						if      (strcmp(szt, "free"      ) == 0) pm->m_bc[bc] =  0;
						else if (strcmp(szt, "fixed"     ) == 0) pm->m_bc[bc] = -1;
						else if (strcmp(szt, "prescribed") == 0)
						{
							pm->m_bc[bc] = lc;
							FERigidBodyDisplacement DC;
							DC.id = nmat;
							DC.bc = bc;
							DC.lc = lc;
							tag.value(DC.sf);
							fem.m_RDC.add(DC);

							// add this boundary condition to the current step
							if (m_nsteps > 0)
							{
								int n = fem.m_RDC.size()-1;
								FERigidBodyDisplacement* pDC = &fem.m_RDC[n];
								m_pStep->AddBoundaryCondition(pDC);
								pDC->Deactivate();
							}
						}
						else if (strcmp(szt, "force") == 0)
						{
							pm->m_bc[bc] = 0;
							FERigidBodyForce FC;
							FC.id = nmat;
							FC.bc = bc;
							FC.lc = lc-1;
							tag.value(FC.sf);
							fem.m_RFC.add(FC);

							// add this boundary condition to the current step
							if (m_nsteps > 0)
							{
								int n = fem.m_RFC.size()-1;
								FERigidBodyForce* pFC = &fem.m_RFC[n];
								m_pStep->AddBoundaryCondition(pFC);
								pFC->Deactivate();
							}
						}
						bfound = true;
					}
				}

				// see if we have processed the tag
				if (bfound == false) throw XMLReader::InvalidTag(tag);
			}

			// get the next tag
			++tag;
		}
		while (!tag.isend());

		// read next tag
		++tag;
	}
	while (!tag.isend());

	// assign material pointers for nested materials
	for (int i=0; i<fem.Materials(); ++i)
	{
		FENestedMaterial* pm = dynamic_cast<FENestedMaterial*>(fem.GetMaterial(i));
		if (pm)
		{
			// get the ID of the base material
			// note that m_nBaseMat is a one-based variable!
			int nbase = pm->m_nBaseMat - 1;

			// make sure the base ID is valid
			if ((nbase < 0) || (nbase >= fem.Materials()))
			{
				fem.m_log.printbox("INPUT ERROR", "Invalid base material ID for material i+1\n");
				throw XMLReader::Error();
				return false;
			}

			// make sure the base material is a valid material (i.e. an elastic material)
			FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(fem.GetMaterial(nbase));
			if (pme == 0)
			{
				fem.m_log.printbox("INPUT ERROR", "Invalid base material for material i+1\n");
				throw XMLReader::Error();
				return false;
			}

			// set the base material pointer
			pm->m_pBase = pme;
		}
	}

	// all done!
	return true;
}

//-----------------------------------------------------------------------------
//! Reads the Nodes section of the FEBio input file

bool FEFEBioImport::ParseNodeSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	int i;
	FEMesh& mesh = fem.m_mesh;

	// first we need to figure out how many nodes there are
	XMLTag t(tag);
	int nodes = 0;
	++t;
	while (!t.isend()) { nodes++; ++t; }

	// create nodes
	mesh.Create(nodes, 0, 0);

	// read nodal coordinates
	++tag;
	for (i=0; i<nodes; ++i)
	{
		FENode& node = fem.m_mesh.Node(i);
		tag.value(node.m_r0);
		node.m_rt = node.m_r0;

		// set rigid body id
		node.m_rid = -1;

		// open displacement dofs
		node.m_ID[0] = 0;
		node.m_ID[1] = 0;
		node.m_ID[2] = 0;

		// open rotational dofs
		node.m_ID[3] = 0;
		node.m_ID[4] = 0;
		node.m_ID[5] = 0;

		// open pressure dof
		node.m_ID[6] = 0;

		// close the rigid rotational dofs
		node.m_ID[7] = -1;
		node.m_ID[8] = -1;
		node.m_ID[9] = -1;

		// open temperature dof
		node.m_ID[10] = 0;

		++tag;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Reads the Element section from the FEBio input file

bool FEFEBioImport::ParseElementSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	int i, j;

	FEMesh& mesh = fem.m_mesh;

	// first we need to figure out how many elements there are
	XMLTag t(tag); ++t;
	int nbel = 0, nsel = 0;
	while (!t.isend())
	{
		if ((t == "hex8") || (t == "penta6") || (t == "tet4")) ++nbel;
		else if ((t == "quad4") || (t == "tri3")) ++nsel;
		else throw XMLReader::InvalidTag(t);

		++t;
	}

	// create elements
	mesh.Create(0, nbel, nsel);

	// read element data
	++tag;
	int n[8];
	int elems = nbel + nsel;
	int nb = 0, ns = 0, nmat;
	FEElement* pe;
	for (i=0; i<elems; ++i)
	{
		if (tag == "hex8")
		{
			FESolidElement& el = mesh.SolidElement(nb++); pe = &el;
			el.SetType(fem.m_nhex8);
			el.m_nID = i+1;
			tag.value(n,el.Nodes());
			for (j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
			nmat = atoi(tag.AttributeValue("mat"))-1;
			el.SetMatID(nmat);
		}
		else if (tag == "penta6")
		{
			FESolidElement& el = mesh.SolidElement(nb++); pe = &el;
			el.SetType(FE_PENTA);
			el.m_nID = i+1;
			tag.value(n,el.Nodes());
			for (j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
			nmat = atoi(tag.AttributeValue("mat"))-1;
			el.SetMatID(nmat);
		}
		else if (tag == "tet4")
		{
			FESolidElement& el = mesh.SolidElement(nb++); pe = &el;
			el.SetType(FE_TET);
			el.m_nID = i+1;
			tag.value(n,el.Nodes());
			for (j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
			nmat = atoi(tag.AttributeValue("mat"))-1;
			el.SetMatID(nmat);
		}
		else if (tag == "quad4")
		{
			FEShellElement& el = mesh.ShellElement(ns++); pe = &el;
			el.SetType(FE_SHELL_QUAD);
			el.m_nID = i+1;
			tag.value(n,el.Nodes());
			for (j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
			nmat = atoi(tag.AttributeValue("mat"))-1;
			el.SetMatID(nmat);

			// thickness are read in the ElementData section
			el.m_h0[0] = el.m_h0[1] = el.m_h0[2] = el.m_h0[3] = 0.0;
		}
		else if (tag == "tri3")
		{
			FEShellElement& el = mesh.ShellElement(ns++); pe = &el;
			el.SetType(FE_SHELL_TRI);
			el.m_nID = i+1;
			tag.value(n,el.Nodes());
			for (j=0; j<el.Nodes(); ++j) el.m_node[j] = n[j]-1;
			nmat = atoi(tag.AttributeValue("mat"))-1;
			el.SetMatID(nmat);

			// thickness are read in the ElementData section
			el.m_h0[0] = el.m_h0[1] = el.m_h0[2] = 0.0;
		}
		else throw XMLReader::InvalidTag(tag);

		// make sure the material number is valid
		if ((nmat<0) || (nmat >= fem.Materials())) throw InvalidMaterial(i+1);

		// go to next tag
		++tag;
	}

	// assign material point data
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);
		FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
		assert(pmat);
		for (int j=0; j<el.GaussPoints(); ++j) el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
	}

	for (i=0; i<mesh.ShellElements(); ++i)
	{
		FEShellElement& el = mesh.ShellElement(i);
		FEMaterial* pmat = fem.GetMaterial(el.GetMatID());
		assert(pmat);
		for (int j=0; j<el.GaussPoints(); ++j) el.SetMaterialPointData(pmat->CreateMaterialPointData(), j);
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Read the NodeData section from the FEBio input file

bool FEFEBioImport::ParseInitialSection(XMLTag& tag)
{
	if (tag.isleaf()) return true;

	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	// make sure we've read the nodes section
	if (mesh.Nodes() == 0) throw XMLReader::InvalidTag(tag);

	for (int i=0; i<mesh.Nodes(); ++i) mesh.Node(i).m_v0 = vec3d(0,0,0);

	// read nodal data
	++tag;
	do
	{
		if (tag == "velocity")
		{
			++tag;
			do
			{
				if (tag == "node")
				{
					int nid = atoi(tag.AttributeValue("id"))-1;
					vec3d v;
					tag.value(v);
					mesh.Node(nid).m_v0 += v;
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Reads the ElementData section from the FEBio input file

bool FEFEBioImport::ParseElementDataSection(XMLTag& tag)
{
	int i;

	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	int nbel = mesh.SolidElements();
	int nsel = mesh.ShellElements();

	//make sure we've read the element section
	int elems = nbel + nsel;
	if (elems == 0) throw XMLReader::InvalidTag(tag);

	// create the pelem array
	vector<FEElement*> pelem;
	pelem.create(nbel + nsel);
	pelem.zero();

	for (i=0; i<nbel; ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);
		assert(pelem[el.m_nID-1] == 0);
		pelem[el.m_nID-1] = &el;
	}

	for (i=0; i<nsel; ++i)
	{
		FEShellElement& el = mesh.ShellElement(i);
		assert(pelem[el.m_nID-1] == 0);
		pelem[el.m_nID-1] = &el;
	}

	// read additional element data
	++tag;
	do
	{
		// make sure this is an "element" tag
		if (tag != "element") throw XMLReader::InvalidTag(tag);

		// get the element number
		const char* szid = tag.AttributeValue("id");
		int n = atoi(szid)-1;

		// make sure the number is valid
		if ((n<0) || (n>=elems)) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

		// get a pointer to the element
		FEElement* pe = pelem[n];

		vec3d a;
		do
		{
			++tag;

			if (tag == "fiber")
			{
				// read the fiber direction
				tag.value(a);

				// normalize fiber
				a.unit();

				// set up a orthonormal coordinate system
				vec3d b(0,1,0);
				if (fabs(fabs(a*b) - 1) < 1e-7) b = vec3d(0,0,1);
				vec3d c = a^b;
				b = c^a;

				// make sure they are unit vectors
				b.unit();
				c.unit();

				FESolidElement* pbe = dynamic_cast<FESolidElement*> (pe);
				FEShellElement* pse = dynamic_cast<FEShellElement*> (pe);
				if (pbe)
				{
					for (int i=0; i<pbe->GaussPoints(); ++i)
					{
						FEElasticMaterialPoint& pt = *pbe->m_State[i]->ExtractData<FEElasticMaterialPoint>();
						mat3d& m = pt.Q;
						m.zero();
						m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
						m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
						m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
					}
				}
				if (pse)
				{
					for (int i=0; i<pse->GaussPoints(); ++i)
					{
						FEElasticMaterialPoint& pt = *pse->m_State[i]->ExtractData<FEElasticMaterialPoint>();
						mat3d& m = pt.Q;
						m.zero();
						m[0][0] = a.x; m[0][1] = b.x; m[0][2] = c.x;
						m[1][0] = a.y; m[1][1] = b.y; m[1][2] = c.y;
						m[2][0] = a.z; m[2][1] = b.z; m[2][2] = c.z;
					}
				}
			}
			else if (tag == "thickness")
			{
				FEShellElement* pse = dynamic_cast<FEShellElement*> (pe);
				if (pse == 0) throw XMLReader::InvalidTag(tag);

				// read shell thickness
				tag.value(pse->m_h0,pse->Nodes());
			}
		}
		while (!tag.isend());

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Reads the Geometry::Groups section of the FEBio input file

bool FEFEBioImport::ParseGroupSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	FEMesh &mesh = fem.m_mesh;

	if (tag.isleaf()) return true;

	++tag;
	do
	{
		if (tag == "node_set")
		{
			// read the node set
			FENodeSet* pns = new FENodeSet(&mesh);

			int nid;
			tag.AttributeValue("id", nid);
			pns->SetID(nid);

			const char* szname = tag.AttributeValue("name", true);
			if (szname) pns->SetName(szname);

			int i, n = 0, n0, n1, nn;
			char* ch;
			const int N = 5096;
			char szval[N], *sz;

			if (strlen(tag.szvalue()) >= N) throw XMLReader::InvalidValue(tag);

			strcpy(szval , tag.szvalue());

			sz = szval;

			int nread;
			do
			{
				ch = strchr(sz, ',');
				if (ch) *ch = 0;
				nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
				switch (nread)
				{
				case 1:
					n1 = n0;
					nn = 1;
					break;
				case 2:
					nn = 1;
					break;
				case 3:
					break;
				default:
					n0 = 0;
					n1 = -1;
					nn = 1;
				}

				for (i=n0; i<=n1; i += nn) ++n;

				if (ch) *ch = ',';
				sz = ch+1;
			}
			while (ch != 0);

			if (n != 0)
			{
				pns->create(n);

				sz = szval;
				n = 0;
				do
				{
					ch = strchr(sz, ',');
					if (ch) *ch = 0;
					nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
					switch (nread)
					{
					case 1:
						n1 = n0;
						nn = 1;
						break;
					case 2:
						nn = 1;
					}

					for (i=n0; i<=n1; i += nn) (*pns)[n++] = i;
					assert(n <= pns->size());

					if (ch) *ch = ',';
					sz = ch+1;
				}
				while (ch != 0);
			}

			// add the nodeset to the mesh
			mesh.AddNodeSet(pns);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseGeometrySection
//  Parses the geometry section from the xml file
//

bool FEFEBioImport::ParseGeometrySection(XMLTag& tag)
{
	++tag;
	do
	{
		if (tag == "Nodes") ParseNodeSection(tag);
		else if (tag == "Elements") ParseElementSection(tag);
		else if (tag == "ElementData") ParseElementDataSection(tag);
		else if (tag == "Groups") ParseGroupSection(tag);
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseBoundarySection
//  Parses the boundary section from the xml file
//

bool FEFEBioImport::ParseBoundarySection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	int n, bc, lc;
	const char* sz;

	// make sure this tag has children
	if (tag.isleaf()) return true;

	++tag;
	do
	{
		if (tag == "fix")
		{
			// make sure this section does not appear in a step section
			if (m_nsteps != 0) throw XMLReader::InvalidTag(tag);

			// Read the fixed nodes
			++tag;
			do
			{
				n = atoi(tag.AttributeValue("id"))-1;
				FENode& node = fem.m_mesh.Node(n);
				sz = tag.AttributeValue("bc");
				if      (strcmp(sz, "x") == 0) node.m_ID[0] = -1;
				else if (strcmp(sz, "y") == 0) node.m_ID[1] = -1;
				else if (strcmp(sz, "z") == 0) node.m_ID[2] = -1;
				else if (strcmp(sz, "xy") == 0) { node.m_ID[0] = node.m_ID[1] = -1; }
				else if (strcmp(sz, "yz") == 0) { node.m_ID[1] = node.m_ID[2] = -1; }
				else if (strcmp(sz, "xz") == 0) { node.m_ID[0] = node.m_ID[2] = -1; }
				else if (strcmp(sz, "xyz") == 0) { node.m_ID[0] = node.m_ID[1] = node.m_ID[2] = -1; }
				else if (strcmp(sz, "p") == 0) { node.m_ID[6] = -1; }
				else if (strcmp(sz, "u") == 0) node.m_ID[3] = -1;
				else if (strcmp(sz, "v") == 0) node.m_ID[4] = -1;
				else if (strcmp(sz, "w") == 0) node.m_ID[5] = -1;
				else if (strcmp(sz, "uv") == 0) { node.m_ID[3] = node.m_ID[4] = -1; }
				else if (strcmp(sz, "vw") == 0) { node.m_ID[4] = node.m_ID[5] = -1; }
				else if (strcmp(sz, "uw") == 0) { node.m_ID[3] = node.m_ID[5] = -1; }
				else if (strcmp(sz, "uvw") == 0) { node.m_ID[3] = node.m_ID[4] = node.m_ID[5] = -1; }
				else if (strcmp(sz, "t") == 0) node.m_ID[10] = -1;
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);
				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "prescribe")
		{
			// count how many prescibed nodes there are
			int ndis = 0;
			XMLTag t(tag); ++t;
			while (!t.isend()) { ndis++; ++t; }

			// allocate prescribed data
			int nsize = fem.m_DC.size();
			fem.m_DC.setsize(nsize + ndis);

			// read the prescribed data
			++tag;
			for (int i=nsize; i<nsize+ndis; ++i)
			{
				n = atoi(tag.AttributeValue("id"))-1;
				sz = tag.AttributeValue("bc");

				if      (strcmp(sz, "x") == 0) bc = 0;
				else if (strcmp(sz, "y") == 0) bc = 1;
				else if (strcmp(sz, "z") == 0) bc = 2;
				else if (strcmp(sz, "p") == 0) bc = 6;	// GAA
				else if (strcmp(sz, "t") == 0) bc = 10; 
	
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

				sz = tag.AttributeValue("lc", true);
				if (sz == 0) lc = 0;
				else lc = atoi(sz);

				fem.m_DC[i].node = n;
				fem.m_DC[i].bc = bc;
				fem.m_DC[i].lc = lc;
				tag.value(fem.m_DC[i].s);

				// add this boundary condition to the current step
				if (m_nsteps > 0)
				{
					m_pStep->AddBoundaryCondition(&fem.m_DC[i]);
					fem.m_DC[i].Deactivate();
				}

				++tag;
			}
		}
		else if (tag == "force")
		{
			// count how many nodal forces there are
			int ncnf = 0;
			XMLTag t(tag); ++t;
			while (!t.isend()) { ncnf++; ++t; }

			// allocate prescribed data
			int nsize = fem.m_FC.size();
			fem.m_FC.setsize(nsize + ncnf);

			// read the prescribed data
			++tag;
			for (int i=nsize; i<ncnf+nsize; ++i)
			{
				n = atoi(tag.AttributeValue("id"))-1;
				sz = tag.AttributeValue("bc");

				if      (strcmp(sz, "x") == 0) bc = 0;
				else if (strcmp(sz, "y") == 0) bc = 1;
				else if (strcmp(sz, "z") == 0) bc = 2;
				else if (strcmp(sz, "p") == 0) bc = 6;	// GAA
				else if (strcmp(sz, "t") == 0) bc = 10;
				else throw XMLReader::InvalidAttributeValue(tag, "bc", sz);

				sz = tag.AttributeValue("lc", true);
				if (sz == 0) lc = 0;
				else lc = atoi(sz);

				fem.m_FC[i].node = n;
				fem.m_FC[i].bc = bc;
				fem.m_FC[i].lc = lc;
				tag.value(fem.m_FC[i].s);

				// add this boundary condition to the current step
				if (m_nsteps > 0)
				{
					m_pStep->AddBoundaryCondition(&fem.m_FC[i]);
					fem.m_FC[i].Deactivate();
				}

				++tag;
			}
		}
		else if (tag == "pressure")
		{
			const char* sz;
			bool blinear = false;
			sz = tag.AttributeValue("type", true);
			if (sz)
			{
				if (strcmp(sz, "linear") == 0) blinear = true;
				else if (strcmp(sz, "nonlinear") == 0) blinear = false;
				else throw XMLReader::InvalidAttributeValue(tag, "type", sz);
			}

			// count how many pressure cards there are
			int npr = 0;
			XMLTag t(tag); ++t;
			while (!t.isend()) { npr++; ++t; }

			// allocate pressure data
			fem.m_PC.create(npr);
			fem.m_psurf->Create(npr);

			// read the pressure data
			++tag;
			int nf[4], N;
			double s;
			for (int i=0; i<npr; ++i)
			{
				FEPressureLoad& pc = fem.m_PC[i];
				FESurfaceElement& el = fem.m_psurf->Element(i);
				pc.blinear = blinear;

				sz = tag.AttributeValue("lc", true);
				if (sz) pc.lc = atoi(sz); else pc.lc = 0;

				const char* sbc = tag.AttributeValue("bc", true);
				if (sbc)
				{
					if (strcmp(sbc, "p") == 0) pc.bc = 0;
					else if (strcmp(sbc, "t") == 0) pc.bc = 1;
					else throw XMLReader::InvalidAttributeValue(tag, "bc", sbc);
				}

				s  = atof(tag.AttributeValue("scale"));
				pc.s[0] = pc.s[1] = pc.s[2] = pc.s[3] = s;

				if (tag == "quad4") el.SetType(FE_QUAD);
				else if (tag == "tri3") el.SetType(FE_TRI);
				else throw XMLReader::InvalidTag(tag);

				N = el.Nodes();
				tag.value(nf, N);
				for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;

				++tag;
			}
		}
		else if (tag == "contact") ParseContactSection(tag);
		else if (tag == "linear_constraint") ParseConstraints(tag);
		else if (tag == "spring")
		{
			FE_DISCRETE_ELEMENT de;
			int n[2];

			// read spring discrete elements
			++tag;
			do
			{
				if (tag == "node")
				{
					tag.value(n, 2);
					de.n1 = n[0]-1;
					de.n2 = n[1]-1;
				}
				else if (tag == "E") tag.value(de.E);
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());

			fem.m_DE.add(de);
		}
		else if (tag == "rigid_body")
		{
			// currently we only allow this to be specified in the multistep feature
			assert(m_nsteps);

			int id;
			tag.AttributeValue("id", id);

			++tag;
			do
			{
				if (tag == "translate")
				{
					const char* szbc = tag.AttributeValue("bc");
					int bc;
					if      (strcmp(szbc, "x") == 0) bc = 0;
					else if (strcmp(szbc, "y") == 0) bc = 1;
					else if (strcmp(szbc, "z") == 0) bc = 2;
					else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

					int lc=-1;
					tag.AttributeValue("lc", lc);

					FERigidBodyDisplacement DC;
					DC.id = id;
					DC.bc = bc;
					DC.lc = lc+1;
					tag.value(DC.sf);
					fem.m_RDC.add(DC);

					// make sure to free the material BC
					FERigid* pm = dynamic_cast<FERigid*>(fem.GetMaterial(id-1));
					assert(pm);
					pm->m_bc[bc] = 0;

					int n = fem.m_RDC.size()-1;
					FERigidBodyDisplacement* pDC = &fem.m_RDC[n];
					pDC->Deactivate();
					fem.m_pStep->AddBoundaryCondition(pDC);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------
//! Parse the linear constraints section of the xml input file
//! This section is a subsection of the Boundary section

bool FEFEBioImport::ParseConstraints(XMLTag& tag)
{
	FEM& fem = *m_pfem;

	// make sure there is a constraint defined
	if (tag.isleaf()) return true;

	// read the master node
	FELinearConstraint LC;
	int node;
	tag.AttributeValue("node", node);
	LC.master.node = node-1;

	const char* szbc = tag.AttributeValue("bc");
	if      (strcmp(szbc, "x") == 0) LC.master.bc = 0;
	else if (strcmp(szbc, "y") == 0) LC.master.bc = 1;
	else if (strcmp(szbc, "z") == 0) LC.master.bc = 2;
	else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

	// we must deactive the master dof
	// so that it does not get assigned an equation
	fem.m_mesh.Node(node-1).m_ID[LC.master.bc] = -1;

	// read the slave nodes
	++tag;
	do
	{
		FELinearConstraint::SlaveDOF dof;
		if (tag == "node")
		{
			tag.value(dof.val);
			tag.AttributeValue("id", node);
			dof.node = node - 1;

			const char* szbc = tag.AttributeValue("bc");
			if      (strcmp(szbc, "x") == 0) dof.bc = 0;
			else if (strcmp(szbc, "y") == 0) dof.bc = 1;
			else if (strcmp(szbc, "z") == 0) dof.bc = 2;
			else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

			LC.slave.push_back(dof);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	// add the linear constraint to the system
	fem.m_LinC.push_back(LC);

	return true;
}


//-----------------------------------------------------------------------------
//! Parses the contact section of the xml input file
//! The contact section is a subsection of the boundary section

bool FEFEBioImport::ParseContactSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	FEMesh& m = fem.m_mesh;

	const char* szt = tag.AttributeValue("type");

	if (strcmp(szt, "sliding_with_gaps") == 0)
	{
		// --- S L I D I N G   W I T H   G A P S ---

		FESlidingInterface* ps = new FESlidingInterface(&fem);
		fem.m_CI.add(ps);

		ps->m_npass = 1;

		++tag;
		do
		{
			if (tag == "tolerance") tag.value(ps->m_atol);
			else if (tag == "gaptol") tag.value(ps->m_gtol);
			else if (tag == "fric_coeff") tag.value(ps->m_mu);
			else if (tag == "fric_penalty") tag.value(ps->m_epsf);
			else if (tag == "minaug") tag.value(ps->m_naugmin);
			else if (tag == "maxaug") tag.value(ps->m_naugmax);
			else if (tag == "search_tol") tag.value(ps->m_stol);
			else if (tag == "ktmult") tag.value(ps->m_ktmult);
			else if (tag == "knmult") tag.value(ps->m_knmult);
			else if (tag == "node_reloc") tag.value(ps->m_breloc);
			else if (tag == "penalty")
			{
				const char* sz = tag.AttributeValue("lc", true);
				if (sz)	ps->m_nplc = atoi(sz);

				sz = tag.AttributeValue("auto", true);
				if (sz)
				{
					if (strcmp(sz, "on") == 0) ps->m_bautopen = true;
				}

				tag.value(ps->m_eps);
			}
			else if (tag == "two_pass")
			{
				int n;
				tag.value(n);
				if ((n<0) || (n>1)) throw XMLReader::InvalidValue(tag);

				ps->m_npass = n+1;
			}
			else if (tag == "seg_up") tag.value(ps->m_nsegup);
			else if (tag == "laugon") tag.value(ps->m_blaugon);
			else if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FEContactSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt);
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "facet-to-facet sliding") == 0)
	{
		// --- F A C E T   T O   F A C E T   S L I D I N G ---

		FEFacet2FacetSliding* ps = new FEFacet2FacetSliding(&fem);
		fem.m_CI.add(ps);

		++tag;
		do
		{
			if (tag == "penalty") tag.value(ps->m_epsn);
			else if (tag == "laugon") tag.value(ps->m_blaugon);
			else if (tag == "tolerance") tag.value(ps->m_atol);
			else if (tag == "gaptol") tag.value(ps->m_gtol);
			else if (tag == "minaug") tag.value(ps->m_naugmin);
			else if (tag == "maxaug") tag.value(ps->m_naugmax);
			else if (tag == "knmult") tag.value(ps->m_knmult);
			else if (tag == "search_tol") tag.value(ps->m_stol);
			else if (tag == "two_pass")
			{
				int n;
				tag.value(n);
				if ((n<0) || (n>1)) throw XMLReader::InvalidValue(tag);

				ps->m_npass = n+1;
			}
			else if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FEFacetSlidingSurface& s = (ntype == 1? ps->m_ms : ps->m_ss);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt);

				// currently the element types are automatically set to FE_NIQUAD or FE_NITRI
				// so we have to modify those elements to FE_QUAD and FE_TRI
				// TODO: we need a better way of doing this!
				for (int i=0; i<s.Elements(); ++i)
				{
					FESurfaceElement& e = s.Element(i);
					if (e.Nodes() == 4) e.SetType(FE_QUAD); 
					else e.SetType(FE_TRI);
				}
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "sliding2") == 0)
	{
		// --- S L I D I N G   I N T E R F A C E   2 ---
		FESlidingInterface2* ps = new FESlidingInterface2(&fem);
		fem.m_CI.add(ps);

		++tag;
		do
		{
			if      (tag == "laugon"    ) tag.value(ps->m_blaugon);
			else if (tag == "tolerance" ) tag.value(ps->m_atol);
			else if (tag == "penalty"   ) tag.value(ps->m_eps);
			else if (tag == "knmult"    ) tag.value(ps->m_knmult);
			else if (tag == "search_tol") tag.value(ps->m_stol);
			else if (tag == "pressure_penalty") tag.value(ps->m_epsp);
			else if (tag == "symmetric_stiffness") tag.value(ps->m_bsymm);
			else if (tag == "search_radius") tag.value(ps->m_srad);
			else if (tag == "two_pass"  ) 
			{
				int n;
				tag.value(n);
				if ((n<0) || (n>1)) throw XMLReader::InvalidValue(tag);

				ps->m_npass = n+1;
			}
			else if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FEContactSurface2& s = (ntype == 1? ps->m_ms : ps->m_ss);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt);

				// currently the element types are automatically set to FE_NIQUAD or FE_NITRI
				// For this type of contact we want gaussian quadrature,
				// so we have to modify those elements to FE_QUAD and FE_TRI
				// TODO: we need a better way of doing this!
				for (int i=0; i<s.Elements(); ++i)
				{
					FESurfaceElement& e = s.Element(i);
					if (e.Nodes() == 4) e.SetType(FE_QUAD); 
					else e.SetType(FE_TRI);
				}
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "tied") == 0)
	{
		// --- T I E D   C O N T A C T  ---

		FETiedInterface* ps = new FETiedInterface(&fem);
		fem.m_CI.add(ps);

		++tag;
		do
		{
			if (tag == "tolerance") tag.value(ps->m_atol);
			else if (tag == "laugon") tag.value(ps->m_blaugon);
			else if (tag == "penalty")
			{
				const char* sz = tag.AttributeValue("lc", true);
				if (sz)	ps->m_nplc = atoi(sz);

				tag.value(ps->m_eps);
			}
			else if (tag == "surface")
			{
				const char* sztype = tag.AttributeValue("type");
				int ntype;
				if (strcmp(sztype, "master") == 0) ntype = 1;
				else if (strcmp(sztype, "slave") == 0) ntype = 2;

				FETiedContactSurface& s = (ntype == 1? ps->ms : ps->ss);

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt);
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "rigid_wall") == 0)
	{
		// --- R I G I D   W A L L   I N T E R F A C E ---

		FERigidWallInterface* ps = new FERigidWallInterface(&fem);
		fem.m_CI.add(ps);

		++tag;
		do
		{
			if (tag == "tolerance") tag.value(ps->m_atol);
			else if (tag == "laugon") tag.value(ps->m_blaugon);
			else if (tag == "penalty")
			{
				const char* sz = tag.AttributeValue("lc", true);
				if (sz)	ps->m_nplc = atoi(sz);

				tag.value(ps->m_eps);
			}
			else if (tag == "plane")
			{
				ps->SetMasterSurface(new FEPlane(&fem));
				FEPlane& pl = dynamic_cast<FEPlane&>(*ps->m_mp);
				const char* sz = tag.AttributeValue("lc", true);
				if (sz)	pl.m_nplc = atoi(sz);

				double* a = pl.GetEquation();
				tag.value(a, 4);
			}
			else if (tag == "sphere")
			{
				ps->SetMasterSurface(new FERigidSphere(&fem));
				FERigidSphere& s = dynamic_cast<FERigidSphere&>(*ps->m_mp);
				++tag;
				do
				{
					if      (tag == "center") tag.value(s.m_rc);
					else if (tag == "radius") tag.value(s.m_R);
					else if (tag == "xtrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[0] = atoi(szlc);
					}
					else if (tag == "ytrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[1] = atoi(szlc);
					}
					else if (tag == "ztrans")
					{
						const char* szlc = tag.AttributeValue("lc");
						s.m_nplc[2] = atoi(szlc);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());
			}
			else if (tag == "surface")
			{
				FEContactSurface& s = ps->m_ss;

				int nfmt = 0;
				const char* szfmt = tag.AttributeValue("format", true);
				if (szfmt)
				{
					if (strcmp(szfmt, "face nodes") == 0) nfmt = 0;
					else if (strcmp(szfmt, "element face") == 0) nfmt = 1;
				}

				// read the surface section
				ParseSurfaceSection(tag, s, nfmt);
			}
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "rigid") == 0)
	{
		// --- R I G I D   B O D Y   I N T E R F A C E ---

		// count how many rigid nodes there are
		int nrn= 0;
		XMLTag t(tag); ++t;
		while (!t.isend()) { nrn++; ++t; }

		// allocate rigid node data
		int nsize = fem.m_RN.size();
		fem.m_RN.setsize(nsize + nrn);

		++tag;
		int id, rb;
		for (int i=nsize; i<nsize+nrn; ++i)
		{
			id = atoi(tag.AttributeValue("id"))-1;
			rb = atoi(tag.AttributeValue("rb"))-1;

			fem.m_RN[i].nid = id;
			fem.m_RN[i].rid = rb;

			if (m_nsteps > 0)
			{
				m_pStep->AddBoundaryCondition(&fem.m_RN[i]);
				fem.m_RN[i].Deactivate();
			}

			++tag;
		}
	}
	else if (strcmp(szt, "rigid joint") == 0)
	{
		// --- R I G I D   J O I N T   I N T E R F A C E ---

		fem.m_nrj++;
		FERigidJoint* prj = new FERigidJoint(&fem);

		++tag;
		do
		{
			if (tag == "tolerance") tag.value(prj->m_atol);
			else if (tag == "penalty") tag.value(prj->m_eps);
			else if (tag =="body_a") tag.value(prj->m_nRBa);
			else if (tag =="body_b") tag.value(prj->m_nRBb);
			else if (tag == "joint") tag.value(prj->m_q0);
			else throw XMLReader::InvalidTag(tag);

			++tag;
		}
		while (!tag.isend());
		prj->m_nRBa--;
		prj->m_nRBb--;
		fem.m_RJ.add(prj);
	}
	else if (strcmp(szt, "linear constraint") == 0)
	{
		FEM& fem = *m_pfem;

		// make sure there is a constraint defined
		if (tag.isleaf()) return true;

		// create a new linear constraint manager
		FELinearConstraintSet* pLCS = new FELinearConstraintSet(&fem);
		fem.m_LCSet.push_back(pLCS);

		// read the linear constraints
		++tag;
		do
		{
			if (tag == "linear_constraint")
			{
				FEAugLagLinearConstraint* pLC = new FEAugLagLinearConstraint;

				FEAugLagLinearConstraint::DOF dof;
				++tag;
				do
				{
					if (tag == "node")
					{
						tag.value(dof.val);
						int node;
						tag.AttributeValue("id", node);
						dof.node = node - 1;

						const char* szbc = tag.AttributeValue("bc");
						if      (strcmp(szbc, "x") == 0) dof.bc = 0;
						else if (strcmp(szbc, "y") == 0) dof.bc = 1;
						else if (strcmp(szbc, "z") == 0) dof.bc = 2;
						else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

						pLC->m_dof.push_back(dof);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());

				// add the linear constraint to the system
				pLCS->add(pLC);
			}
			else if (tag == "tol"    ) tag.value(pLCS->m_tol);
			else if (tag == "penalty") tag.value(pLCS->m_eps);
			else if (tag == "maxaug") tag.value(pLCS->m_naugmax);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else if (strcmp(szt, "traction constraint") == 0)
	{
		FEM& fem = *m_pfem;

		// make sure there is a constraint defined
		if (tag.isleaf()) return true;

		// create a new traction constraint manager
		FETractionConstraintSet* pRCS = new FETractionConstraintSet(&fem);
		fem.m_RCSet.push_back(pRCS);

		// read the traction constraints
		++tag;
		do
		{
			if (tag == "traction_constraint")
			{
				FETractionConstraint* pRC = new FETractionConstraint;

				FETractionConstraint::DOF dof;
				++tag;
				do
				{
					if (tag == "node")
					{
						tag.value(dof.val);
						int node;
						tag.AttributeValue("id", node);
						dof.node = node - 1;

						const char* szbc = tag.AttributeValue("bc");
						if      (strcmp(szbc, "x") == 0) dof.bc = 0;
						else if (strcmp(szbc, "y") == 0) dof.bc = 1;
						else if (strcmp(szbc, "z") == 0) dof.bc = 2;
						else throw XMLReader::InvalidAttributeValue(tag, "bc", szbc);

						pRC->m_dof.push_back(dof);
					}
					else throw XMLReader::InvalidTag(tag);
					++tag;
				}
				while (!tag.isend());

				// add the residual constraint to the system
				pRCS->add(pRC);
			}
			else if (tag == "tol"    ) tag.value(pRCS->m_tol);
			else if (tag == "penalty") tag.value(pRCS->m_eps);
			else if (tag == "maxaug") tag.value(pRCS->m_naugmax);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
		while (!tag.isend());
	}
	else throw XMLReader::InvalidAttributeValue(tag, "type", szt);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseGlobalsSection
//  This function reads the global variables from the xml file
//

bool FEFEBioImport::ParseGlobalsSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;

	fem.m_BF[0].lc = -1;
	fem.m_BF[1].lc = -1;
	fem.m_BF[2].lc = -1;

	++tag;
	do
	{
		if (tag == "body_force")
		{
			++tag;
			int n;
			const char* szlc;
			do
			{
				n = -1;
				if (tag == "x") n = 0;
				else if (tag == "y") n = 1;
				else if (tag == "z") n = 2;

				if (n == -1) throw XMLReader::InvalidTag(tag);

				szlc = tag.AttributeValue("lc");
				fem.m_BF[n].lc = atoi(szlc);
				tag.value(fem.m_BF[n].s);

				++tag;
			}
			while (!tag.isend());
		}

		++tag;
	}
	while (!tag.isend());

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION: FEFEBioImport::ParseLoadSection
//  This function reads the load data section from the xml file
//

bool FEFEBioImport::ParseLoadSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;

	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			// count how many points we have
			XMLTag t(tag); ++t;
			int nlp = 0;
			while (!t.isend()) { ++nlp; ++t; }

			// create the loadcurve
			FELoadCurve* plc = new FELoadCurve;
			plc->Create(nlp);
			fem.AddLoadCurve(plc);

			// read the points
			double d[2];
			++tag;
			for (int i=0; i<nlp; ++i)
			{
				tag.value(d, 2);
				plc->LoadPoint(i).time  = d[0];
				plc->LoadPoint(i).value = d[1];

				++tag;
			}
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//-----------------------------------------------------------------------------

bool FEFEBioImport::ParseOutputSection(XMLTag& tag)
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	const char* sz;

	++tag;
	do
	{
		if (tag == "logfile")
		{
			++tag;
			do
			{
				if (tag == "node_data")
				{
					sz = tag.AttributeValue("file", true);

					NodeDataRecord* prec = new NodeDataRecord(m_pfem, sz);
					sz = tag.AttributeValue("data");
					strcpy(prec->m_szdata, sz);

					sz = tag.AttributeValue("name", true);
					if (sz != 0) strcpy(prec->m_szname, sz);

					sz = tag.AttributeValue("delim", true);
					if (sz != 0) strcpy(prec->m_szdelim, sz);

					sz = tag.AttributeValue("comments", true);
					if (sz != 0)
					{
						if      (strcmp(sz, "on") == 0) prec->m_bcomm = true;
						else if (strcmp(sz, "off") == 0) prec->m_bcomm = false; 
					}

					if (tag.isleaf()) prec->DataRecord::SetItemList(tag.szvalue());
					else
					{
						++tag;
						if (tag == "node_set")
						{
							FENodeSet* pns = 0;
							const char* szid = tag.AttributeValue("id", true);
							if (szid == 0)
							{
								const char* szname = tag.AttributeValue("name");
								pns = mesh.FindNodeSet(szname);
							}
							else pns = mesh.FindNodeSet(atoi(szid));

							if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, "id", szid);

							prec->SetItemList(pns);
						}
						else throw XMLReader::InvalidTag(tag);
						++tag;
						assert(tag.isend());
					}

					fem.m_Data.AddRecord(prec);
				}
				else if (tag == "element_data")
				{
					sz = tag.AttributeValue("file", true);

					ElementDataRecord* prec = new ElementDataRecord(m_pfem, sz);
					sz = tag.AttributeValue("data");
					strcpy(prec->m_szdata, sz);

					sz = tag.AttributeValue("name", true);
					if (sz != 0) strcpy(prec->m_szname, sz);


					sz = tag.AttributeValue("delim", true);
					if (sz != 0) strcpy(prec->m_szdelim, sz);

					sz = tag.AttributeValue("comments", true);
					if (sz != 0)
					{
						if      (strcmp(sz, "on") == 0) prec->m_bcomm = true;
						else if (strcmp(sz, "off") == 0) prec->m_bcomm = false; 
					}

					prec->SetItemList(tag.szvalue());

					fem.m_Data.AddRecord(prec);
				}
				else if (tag == "rigid_body_data")
				{
					sz = tag.AttributeValue("file", true);

					RigidBodyDataRecord* prec = new RigidBodyDataRecord(m_pfem, sz);
					sz = tag.AttributeValue("data");
					strcpy(prec->m_szdata, sz);

					sz = tag.AttributeValue("name", true);
					if (sz != 0) strcpy(prec->m_szname, sz);

					sz = tag.AttributeValue("delim", true);
					if (sz != 0) strcpy(prec->m_szdelim, sz);

					sz = tag.AttributeValue("comments", true);
					if (sz != 0)
					{
						if      (strcmp(sz, "on") == 0) prec->m_bcomm = true;
						else if (strcmp(sz, "off") == 0) prec->m_bcomm = false; 
					}

					prec->SetItemList(tag.szvalue());

					fem.m_Data.AddRecord(prec);
				}
				else throw XMLReader::InvalidTag(tag);

				++tag;
			}
			while (!tag.isend());
		}
		else if (tag == "plotfile")
		{
			++tag;
			do
			{
				if (tag == "shell_strain") tag.value(fem.m_plot.m_bsstrn);
				else if (tag == "map")
				{
					const char* szfield = tag.AttributeValue("field");
					const char* szval = tag.szvalue();
					if (strcmp(szfield, "displacement") == 0)
					{
						if (strcmp(szval, "DISPLACEMENT") == 0) fem.m_plot.m_nfield[0] = PlotFile::PLOT_DISPLACEMENT;
						else throw XMLReader::InvalidValue(tag);
					}
					else if (strcmp(szfield, "velocity") == 0)
					{
						if (strcmp(szval, "NONE") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_NONE;
						else if (strcmp(szval, "VELOCITY") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_VELOCITY;
						else if (strcmp(szval, "FLUID_FLUX") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_FLUID_FLUX;
						else if (strcmp(szval, "CONTACT_TRACTION") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_CONTACT_TRACTION;
						else if (strcmp(szval, "REACTION_FORCE") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_REACTION_FORCE;
						else if (strcmp(szval, "MATERIAL_FIBER") == 0) fem.m_plot.m_nfield[1] = PlotFile::PLOT_MATERIAL_FIBER;
						else throw XMLReader::InvalidValue(tag);
					}
					else if (strcmp(szfield, "acceleration") == 0)
					{
						if (strcmp(szval, "NONE") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_NONE;
						else if (strcmp(szval, "ACCELERATION") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_ACCELERATION;
						else if (strcmp(szval, "FLUID_FLUX") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_FLUID_FLUX;
						else if (strcmp(szval, "CONTACT_TRACTION") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_CONTACT_TRACTION;
						else if (strcmp(szval, "REACTION_FORCE") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_REACTION_FORCE;
						else if (strcmp(szval, "MATERIAL_FIBER") == 0) fem.m_plot.m_nfield[2] = PlotFile::PLOT_MATERIAL_FIBER;
						else throw XMLReader::InvalidValue(tag);
					}
					else if (strcmp(szfield, "temperature") == 0)
					{
						if (strcmp(szval, "NONE") == 0) fem.m_plot.m_nfield[3] = PlotFile::PLOT_NONE;
						else if (strcmp(szval, "FLUID_PRESSURE") == 0) fem.m_plot.m_nfield[3] = PlotFile::PLOT_FLUID_PRESSURE;
						else if (strcmp(szval, "CONTACT_PRESSURE") == 0) fem.m_plot.m_nfield[3] = PlotFile::PLOT_CONTACT_PRESSURE;
						else if (strcmp(szval, "CONTACT_GAP") == 0) fem.m_plot.m_nfield[3] = PlotFile::PLOT_CONTACT_GAP;
						else throw XMLReader::InvalidValue(tag);
					}
					else if (strcmp(szfield, "plastic strain") == 0)
					{
						if      (strcmp(szval, "PLASTIC_STRAIN"  ) == 0) fem.m_plot.m_nfield[4] = PlotFile::PLOT_PLASTIC_STRAIN;
						else if (strcmp(szval, "FIBER_STRAIN"    ) == 0) fem.m_plot.m_nfield[4] = PlotFile::PLOT_FIBER_STRAIN;
						else if (strcmp(szval, "DEV_FIBER_STRAIN") == 0) fem.m_plot.m_nfield[4] = PlotFile::PLOT_DEV_FIBER_STRAIN;
						else throw XMLReader::InvalidValue(tag);
					}
					else throw XMLReader::InvalidAttributeValue(tag, "field", szfield);
				}
				else throw XMLReader::InvalidTag(tag);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());

	return  true;
}

//-----------------------------------------------------------------------------

bool FEFEBioImport::ParseStepSection(XMLTag& tag)
{
	// We assume that the FEM object will already have at least one step
	// defined. Therefor the first time we find a "step" section we
	// do not create a new step. If more steps are required
	// we need to create new FEAnalysis steps and add them to fem
	if (m_nsteps != 0)
	{
		m_pStep = new FEAnalysis(*m_pfem);
		m_pfem->m_Step.add(m_pStep);
	}

	// increase the step section counter
	++m_nsteps;

	++tag;
	do
	{
		if (tag == "Control") ParseControlSection(tag);
		else if (tag == "Boundary") ParseBoundarySection(tag);
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());

	return true;
}

//---------------------------------------------------------------------------------

bool FEFEBioImport::ParseSurfaceSection(XMLTag &tag, FESurface& s, int nfmt)
{
	FEM& fem = *m_pfem;
	FEMesh& m = fem.m_mesh;

	// count nr of faces
	int faces = 0, N, nf[4];
	XMLTag t(tag); ++t;
	while (!t.isend()) { faces++; ++t; }

	// allocate storage for faces
	s.Create(faces);

	// read faces
	++tag;
	for (int i=0; i<faces; ++i)
	{
		FESurfaceElement& el = s.Element(i);

		if (tag == "quad4") el.SetType(FE_NIQUAD);
		else if (tag == "tri3") el.SetType(FE_NITRI);
		else throw XMLReader::InvalidTag(tag);

		N = el.Nodes();

		if (nfmt == 0)
		{
			tag.value(nf, N);
			for (int j=0; j<N; ++j) el.m_node[j] = nf[j]-1;
		}
		else if (nfmt == 1)
		{
			tag.value(nf, 2);
			FEElement* pe = m.FindElementFromID(nf[0]);
			if (pe)
			{
				int ne[4];
				int nn = m.GetFace(*pe, nf[1]-1, ne);
				if (nn != N) throw XMLReader::InvalidValue(tag);
				for (int j=0; j<N; ++j) el.m_node[j] = ne[j];
			}
			else throw XMLReader::InvalidValue(tag);
		}

		++tag;
	}
	return true;
}
