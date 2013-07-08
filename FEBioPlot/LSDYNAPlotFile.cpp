#include "stdafx.h"
#include "LSDYNAPlotFile.h"
#include "FECore/FEDiscreteDomain.h"
#include "FECore/FETrussDomain.h"
#include "FECore/FEShellDomain.h"
#include "FECore/FESolidDomain.h"
#include "FEBioLib/FESolidSolver.h"
#include "FEBioLib/FESlidingInterface.h"
#include "FEBioLib/FETiedInterface.h"
#include "FEBioLib/FERigidWallInterface.h"
#include "FEBioLib/FEPeriodicBoundary.h"
#include "FEBioLib/FESurfaceConstraint.h"
#include "FEBioLib/FESlidingInterface2.h"
#include "FEBioLib/FESlidingInterface3.h"
#include "FEBioLib/FEFacet2FacetSliding.h"
#include "FEBioLib/FETransverselyIsotropic.h"
#include "FEBioLib/FEHeatTransferMaterial.h"
#include "FEBioLib/FETrussMaterial.h"
#include "FEBioLib/FEBiphasic.h"
#include "FEBioLib/FEAnalysisStep.h"

LSDYNAPlotFile::LSDYNAPlotFile()
{
	m_bsstrn = false;
	m_nfield[0] = -1;	// displacement
	m_nfield[1] = -1;	// velocity
	m_nfield[2] = -1;	// acceleration
	m_nfield[3] = -1;	// temperature
	m_nfield[4] = -1;	// plastic strain

	m_fp = 0;
}

LSDYNAPlotFile::~LSDYNAPlotFile()
{
	if (m_fp) fclose(m_fp);
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : PlotFile::Open
//  Opens the PLOT database. That is, an archive is created and the PLOTHEADER
//  is written as well as the initial geometry of the model.
//

bool LSDYNAPlotFile::Open(FEModel& fem, const char* szfile)
{
	int i, j, N, nd;

	// open the archive
	if ((m_fp = fopen(szfile, "wb")) == 0) return false;

	// store a pointer to the FEModel class
	m_pfem = &fem;

	// get the mesh
	FEMesh& mesh = fem.GetMesh();

	int nmode = fem.GetCurrentStep()->GetType();
	int ntype = fem.GetCurrentStep()->m_nanalysis;

	// check the field values
	if (m_nfield[0] == -1) 
	{
//		if (itype == FE_HEAT_CONDUCTION) m_nfield[0] = PLOT_NONE;
//		else m_nfield[0] = PLOT_DISPLACEMENT;
		m_nfield[0] = PLOT_DISPLACEMENT;
	}
	if (m_nfield[1] == -1)
	{
		m_nfield[1] = PLOT_NONE;
		if ((nmode == FE_BIPHASIC) || (ntype == FE_DYNAMIC)) m_nfield[1] = PLOT_VELOCITY;
		else if (nmode == FE_HEAT) m_nfield[1] = PLOT_HEAT_FLUX;
	}
	if (m_nfield[2] == -1)
	{
		m_nfield[2] = PLOT_NONE;
		if (nmode == FE_BIPHASIC) m_nfield[2] = PLOT_FLUID_FLUX;
		else if (ntype == FE_DYNAMIC) m_nfield[2] = PLOT_ACCELERATION;
		else if (fem.ContactInterfaces() > 0) m_nfield[2] = PLOT_CONTACT_TRACTION;
	}
	if (m_nfield[3] == -1)
	{
		m_nfield[3] = PLOT_NONE;
		if (nmode == FE_BIPHASIC) m_nfield[3] = PLOT_FLUID_PRESSURE;
		else if (nmode == FE_HEAT) m_nfield[3] = PLOT_TEMPERATURE;
		else if (fem.ContactInterfaces() > 0) m_nfield[3] = PLOT_CONTACT_GAP;
	}
	if (m_nfield[4] == -1)
	{
		m_nfield[4] = PLOT_DEV_FIBER_STRAIN;
	}

	// write the header
	PLOTHEADER plh = {0};

	// copy the title
	// note that we only have 40 characters to store the title
	const char* sztitle = fem.GetTitle();
	if (strlen(sztitle) > 39)
		strncpy(plh.Title, sztitle, 39);
	else
		strcpy(plh.Title, fem.GetTitle());

	plh.neips = 2000;
	plh.flagU = (m_nfield[0] == 0? 0 : 1);
	plh.flagV = (m_nfield[1] == 0? 0 : 1);
	plh.flagA = (m_nfield[2] == 0? 0 : 1);
	plh.flagT = (m_nfield[3] == 0? 0 : 1);
	plh.icode = 6;
	plh.ndim  = 4;
	plh.nel2  = mesh.TrussElements() + mesh.DiscreteElements();
	plh.nel4  = mesh.ShellElements();
	plh.nel8  = mesh.SolidElements();
	plh.nglbv = 0;
	plh.nummat2 = (mesh.DiscreteElements()>0? 1 : 0);
	plh.nummat4 = 0;
	plh.nummat8 = fem.Materials();
	plh.nump    = mesh.Nodes();
	plh.nv1d    = 6;
	plh.nv2d    = (m_bsstrn? 44 : 32);
	plh.nv3d    = 7;
	plh.maxint  = 0;
	plh.neips   = 2000;
	plh.ioshl4  = 0;

	// store the header
	m_ph = plh;

	fwrite(&plh, sizeof(PLOTHEADER), 1, m_fp);

	// write the material coordinates
	float xf[3];
	for (i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		xf[0] = (float) node.m_r0.x;
		xf[1] = (float) node.m_r0.y;
		xf[2] = (float) node.m_r0.z;

		fwrite(xf, sizeof(float), 3, m_fp);
	}

	// write the connectivity and material number
	// Note that we increment all numbers by 1 since
	// the plot database expects 1-based arrays
	int n[9];

	int nid = 1;

	// write solid element data
	// note that we reindex all elements so that the ID
	// corresponds to the nr in the plot file
	for (nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				el.m_nID = nid++;

				N = el.Nodes();
				switch (el.Type())
				{
				case FE_HEX8G8:
				case FE_HEX8RI:
				case FE_HEX8G1:
					for (j=0; j<N; ++j) n[j] = el.m_node[j]+1;
					break;
				case FE_PENTA6G6:
					// note the weird mapping. This is to be consistent
					// with NIKE's wedge element
					n[0] = el.m_node[0]+1;
					n[1] = el.m_node[2]+1;
					n[2] = el.m_node[5]+1;
					n[3] = el.m_node[3]+1;
					n[4] = el.m_node[1]+1;
					n[5] = el.m_node[1]+1;
					n[6] = el.m_node[4]+1;
					n[7] = el.m_node[4]+1;
					break;
				case FE_TET4G4:
				case FE_TET4G1:
					n[0] = el.m_node[0]+1;
					n[1] = el.m_node[1]+1;
					n[2] = el.m_node[2]+1;
					n[3] = el.m_node[2]+1;
					n[4] = n[5] = n[6] = n[7] = el.m_node[3]+1;
					break;
				}

				n[8] = el.GetMatID()+1;

				fwrite(n, sizeof(int), 9, m_fp);
			}
		}
	}

	// write truss element data
	for (nd=0; nd < mesh.Domains(); ++nd)
	{
		FETrussDomain* ptd = dynamic_cast<FETrussDomain*>(&mesh.Domain(nd));
		if (ptd)
		{
			for (i=0; i<ptd->Elements(); ++i)
			{
				FETrussElement& el = ptd->Element(i);
				el.m_nID = nid++;
				n[0] = el.m_node[0]+1;
				n[1] = el.m_node[1]+1;
				n[2] = 0;
				n[3] = 0;
				n[4] = 0;
				n[5] = el.GetMatID()+1;

				fwrite(n, sizeof(int), 6, m_fp);
			}
		}
	}

	// all springs are assigned the same material for now
	int nmat = fem.Materials()+1;
	for (nd=0; nd < mesh.Domains(); ++nd)
	{
		FEDiscreteDomain* psd = dynamic_cast<FEDiscreteDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEDiscreteElement& e = dynamic_cast<FEDiscreteElement&>(psd->ElementRef(i));
				n[0] = e.m_node[0]+1;
				n[1] = e.m_node[1]+1;
				n[2] = 0;
				n[3] = 0;
				n[4] = 0;
				n[5] = nmat;

				fwrite(n, sizeof(int), 6, m_fp);
			}
		}
	}

	// write shell element data
	for (nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);

				el.m_nID = nid++;

				N = el.Nodes();
				switch (el.Type())
				{
				case FE_SHELL_QUAD:
					n[0] = el.m_node[0]+1;
					n[1] = el.m_node[1]+1;
					n[2] = el.m_node[2]+1;
					n[3] = el.m_node[3]+1;
					break;
				case FE_SHELL_TRI:
					n[0] = el.m_node[0]+1;
					n[1] = el.m_node[1]+1;
					n[2] = el.m_node[2]+1;
					n[3] = el.m_node[2]+1;
					break;
				}
		
				n[4] = el.GetMatID()+1;
		
				fwrite(n, sizeof(int), 5, m_fp);
			}
		}
	}

	return true;
}


//-----------------------------------------------------------------------------
//! Opens a PLOT database for appending
//!
//! \todo I always need to make sure that the same amount of data is written
//!       to the file after opening for appending. In other words, make sure that
//!       the flags are the same as when opening the file for the first time.
bool LSDYNAPlotFile::Append(FEModel& fem, const char* szfile)
{
	// open the file
	if ((m_fp = fopen(szfile, "rb")) == 0) return false;
	
	// read the plot file header
	fread(&m_ph, sizeof(PLOTHEADER), 1, m_fp);

	if (m_ph.nv2d == 32) m_bsstrn = false;
	else if (m_ph.nv2d == 44) m_bsstrn = true;
	else
	{
		fclose(m_fp); m_fp = 0;
		return false;
	}

	// close the file
	fclose(m_fp);

	m_pfem = &fem;

	// reopen the plot file for appending
	m_fp = fopen(szfile, "a+b");
	return (m_fp != 0);
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : PlotFile::Write
//  Writes the current time step data to the PLOT database. The data archive
//  must have been opened with the Open or Append function.
//

bool LSDYNAPlotFile::Write(FEModel& fem)
{
	// make sure the archive is opened
	if (m_fp == 0) return false;

	// write the time value
	float time = (float) fem.m_ftime;
	fwrite(&time, sizeof(float) ,1, m_fp);

	FEMesh& mesh = fem.GetMesh();

	m_pfem = &fem;

	// write the spatial coordinates
	if (m_ph.flagU)
	{
		switch (m_nfield[0])
		{
		case PLOT_DISPLACEMENT: write_displacements(); break;
		default:
			assert(false);
		}
	}

	// write the velocities
	if (m_ph.flagV)
	{
		switch (m_nfield[1])
		{
		case PLOT_VELOCITY: write_velocities(); break;
		case PLOT_FLUID_FLUX: write_fluid_flux(); break;
		case PLOT_CONTACT_TRACTION: write_contact_tractions(); break;
		case PLOT_REACTION_FORCE: write_reaction_forces(); break;
		case PLOT_MATERIAL_FIBER: write_material_fibers(); break;
		case PLOT_HEAT_FLUX: write_heat_flux(); break;
		default:
			assert(false);
		}
	}

	// write the accelerations
	if (m_ph.flagA)
	{
		switch (m_nfield[2])
		{
		case PLOT_ACCELERATION: write_accelerations(); break;
		case PLOT_FLUID_FLUX: write_fluid_flux(); break;
		case PLOT_CONTACT_TRACTION: write_contact_tractions(); break;
		case PLOT_REACTION_FORCE: write_reaction_forces(); break;
		case PLOT_MATERIAL_FIBER: write_material_fibers(); break;
		case PLOT_HEAT_FLUX: write_heat_flux(); break;
		default:
			assert(false);
		}
	}

	// write the temperatures
	if (m_ph.flagT)
	{
		switch (m_nfield[3])
		{
		case PLOT_FLUID_PRESSURE: write_fluid_pressures(); break;
		case PLOT_CONTACT_PRESSURE: write_contact_pressures(); break;
		case PLOT_CONTACT_GAP: write_contact_gaps(); break;
		case PLOT_TEMPERATURE: write_temperatures(); break;
		default:
			assert(false);
		}
	}

	// write the stresses
	// note that at this point the plastic strain field (= s[6]) is not written to
	// or better, is always set to zero.

	// write solid element data
	write_solid_stress();

	// write truss element data
	write_truss_stress();

	// write shell data
	write_shell_stress();

	// now flush the archive to make sure we don't loose any results
	fflush(m_fp);

	return true;
}

//-----------------------------------------------------------------------------
void LSDYNAPlotFile::write_solid_stress()
{
	int i, j;
	float s[7] = {0};
	double f;
	int nint, nd;
	FEMesh& mesh = m_pfem->GetMesh();
	for (nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				for (j=0; j<7; ++j) s[j] = 0;

				nint = el.GaussPoints();

				f = 1.0 / (double) nint;
		
				// since the PLOT file requires floats we need to convert
				// the doubles to single precision
				// we output the average stress values of the gauss points
				for (j=0; j<nint; ++j)
				{
					FEElasticMaterialPoint* ppt = (el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

					if (ppt)
					{
						FEElasticMaterialPoint& pt = *ppt;
						s[0] += (float) (f*pt.s.xx());
						s[1] += (float) (f*pt.s.yy());
						s[2] += (float) (f*pt.s.zz());
						s[3] += (float) (f*pt.s.xy());
						s[4] += (float) (f*pt.s.yz());
						s[5] += (float) (f*pt.s.xz());

						if (m_nfield[4] == PLOT_FIBER_STRAIN)
						{
							s[6] += (float) f*fiber_strain(el, j);
						}
						else if (m_nfield[4] == PLOT_DEV_FIBER_STRAIN)
						{
							s[6] += (float) f*dev_fiber_strain(el, j);
						}
					}
				}

				fwrite(s, sizeof(float), 7, m_fp);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void LSDYNAPlotFile::write_shell_stress()
{
	int i, j, neln;
	float s[44] = {0};

	// write shell element data
	mat3ds E;
	FEMesh& mesh = m_pfem->GetMesh();
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				for (j=0; j<44; ++j) s[j] = 0;

				FEShellElement& el = psd->Element(i);
				if (!el.IsRigid())
				{
					try
					{
						neln = el.Nodes();
						double f = 1.0 / (double) neln;

						// output shell data
						for (j=0; j<neln; ++j)
						{
							FEElasticMaterialPoint& ptm = *(el.m_State[j + neln]->ExtractData<FEElasticMaterialPoint>());
							FEElasticMaterialPoint& pti = *(el.m_State[j       ]->ExtractData<FEElasticMaterialPoint>());
							FEElasticMaterialPoint& pto = *(el.m_State[j+2*neln]->ExtractData<FEElasticMaterialPoint>());

							// mid-surface stresses
							s[ 0] += (float) (f*ptm.s.xx());
							s[ 1] += (float) (f*ptm.s.yy());
							s[ 2] += (float) (f*ptm.s.zz());
							s[ 3] += (float) (f*ptm.s.xy());
							s[ 4] += (float) (f*ptm.s.yz());
							s[ 5] += (float) (f*ptm.s.xz());

							// inner surface stresses
							s[ 7] += (float) (f*pti.s.xx());
							s[ 8] += (float) (f*pti.s.yy());
							s[ 9] += (float) (f*pti.s.zz());
							s[10] += (float) (f*pti.s.xy());
							s[11] += (float) (f*pti.s.yz());
							s[12] += (float) (f*pti.s.xz());

							// outer surface stresses
							s[14] += (float) (f*pto.s.xx());
							s[15] += (float) (f*pto.s.yy());
							s[16] += (float) (f*pto.s.zz());
							s[17] += (float) (f*pto.s.xy());
							s[18] += (float) (f*pto.s.yz());
							s[19] += (float) (f*pto.s.xz());

							// shell thicknesses
							s[29] += (float) (el.m_h0[j]*f*mesh.Node(el.m_node[j]).m_Dt.norm());

							if (m_bsstrn && (!el.IsRigid()))
							{
								// inner-surface strain
								E = pti.Strain();

								s[32] += (float) (f*E.xx());
								s[33] += (float) (f*E.yy());
								s[34] += (float) (f*E.zz());
								s[35] += (float) (f*E.xy());
								s[36] += (float) (f*E.yz());
								s[37] += (float) (f*E.xz());

								// outer-surface strain
								E = pto.Strain();

								s[38] += (float) (f*E.xx());
								s[39] += (float) (f*E.yy());
								s[40] += (float) (f*E.zz());
								s[41] += (float) (f*E.xy());
								s[42] += (float) (f*E.yz());
								s[43] += (float) (f*E.xz());
							}
						}
					}
					catch (...)
					{
						// don't do anything
					}
				} // if (!el.isrigid())

				// save data to file
				fwrite(s, sizeof(float), (m_bsstrn?44:32), m_fp);
			}
		}
	}
}

//----------------------------------------------------------------------------
void LSDYNAPlotFile::write_truss_stress()
{
	float s[6] = {0};
	FEMesh& mesh = m_pfem->GetMesh();
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FETrussDomain* ptd = dynamic_cast<FETrussDomain*>(&mesh.Domain(nd));
		if (ptd)
		{
			for (int i=0; i<ptd->Elements(); ++i)
			{
				FETrussElement& el = ptd->Element(i);
				FETrussMaterialPoint& pt = *(el.m_State[0]->ExtractData<FETrussMaterialPoint>());

				vec3d r0[2], rt[2];
				for (int j=0; j<2; ++j)
				{
					r0[j] = mesh.Node(el.m_node[j]).m_r0;
					rt[j] = mesh.Node(el.m_node[j]).m_rt;
				}

				double l = (rt[1] - rt[0]).norm();
				double L = (r0[1] - r0[0]).norm();
				double V = L*el.m_a0;

				s[0] = (float) (pt.m_tau*V/l);	// axial force

				fwrite(s, sizeof(float), 6, m_fp);
			}
		}
	}

	s[0] = s[1] = s[2] = s[3] = s[4] = s[5] = 0;
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FEDiscreteDomain* psd = dynamic_cast<FEDiscreteDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			for (int i=0; i<psd->Elements(); ++i)
			{
				FEDiscreteElement& e = dynamic_cast<FEDiscreteElement&>(psd->ElementRef(i));
				fwrite(s, sizeof(float), 6, m_fp);
			}
		}
	}
}

//-----------------------------------------------------------------------------

void LSDYNAPlotFile::write_velocities()
{
	FEMesh& mesh = m_pfem->GetMesh();
	float vf[3];
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		vf[0] = (float) node.m_vt.x;
		vf[1] = (float) node.m_vt.y;
		vf[2] = (float) node.m_vt.z;

		fwrite(vf, sizeof(float), 3, m_fp);
	}
}

void LSDYNAPlotFile::write_accelerations()
{
	FEMesh& mesh =m_pfem->GetMesh();
	float af[3];
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		af[0] = (float) node.m_at.x;
		af[1] = (float) node.m_at.y;
		af[2] = (float) node.m_at.z;

		fwrite(af, sizeof(float), 3, m_fp);
	}
}

void LSDYNAPlotFile::write_fluid_flux()
{
	FEMesh& mesh = m_pfem->GetMesh();

	int i, j;

	vector<vec3d> wn(mesh.Nodes());
	vector<int> val(mesh.Nodes());
	for (i=0; i<mesh.Nodes(); ++i) val[i] = 0;

	vec3d ew;
	int n;
	for (int nd=0; nd<mesh.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				// calculate average flux
				ew = vec3d(0,0,0);
				for (j=0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.m_State[j];
					FEBiphasicMaterialPoint* pt = (mp.ExtractData<FEBiphasicMaterialPoint>());

					if (pt) ew += pt->m_w;
				}

				ew /= el.GaussPoints();

				// project to nodes
				for (j=0; j<el.Nodes(); ++j)
				{
					n = el.m_node[j];
					wn[n] += ew;
					++val[n];
				}
			}
		}
	}
	for (i=0; i<mesh.Nodes(); ++i) if (val[i] != 0) wn[i] /= val[i];

	// output nodal fluxes
	float af[3];
	for (i=0; i<mesh.Nodes(); ++i)
	{
		af[0] = (float) wn[i].x;
		af[1] = (float) wn[i].y;
		af[2] = (float) wn[i].z;

		fwrite(af, sizeof(float), 3, m_fp);
	}
}

void LSDYNAPlotFile::write_contact_tractions()
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	int i, j, k, n;

	vector<float> acc(3*mesh.Nodes()); zero(acc);
	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FEContactInterface* pci = fem.ContactInterface(i);

		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*> (pci);
		if (psi)
		{
			int npass = (psi->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FESlidingSurface& ss = (n==0?psi->m_ss:psi->m_ms);
				FESlidingSurface& ms = (n==0?psi->m_ms:psi->m_ss);
				for (j=0; j<ss.Nodes(); ++j)
				{
					int m = ss.m_node[j];
					vec3d t = ss.traction(j);

					acc[3*m  ] += (float) t.x;
					acc[3*m+1] += (float) t.y;
					acc[3*m+2] += (float) t.z;
						
//					acc[m][0] = (float) ss.Lt[j][0];
//					acc[m][1] = (float) ss.Lt[j][1];
//					acc[m][2] = (float) ss.Lm[j];
				}
			}
		}

		FEPeriodicBoundary* pbi = dynamic_cast<FEPeriodicBoundary*>(pci);
		if (pbi)
		{
			FEPeriodicSurface& ss = pbi->m_ss;
			FEPeriodicSurface& ms = pbi->m_ms;
			for (j=0; j<ss.Nodes(); ++j)
			{
				vec3d t = ss.m_Lm[j];// + ss.m_gap[j]*pbi->m_eps;
				int m = ss.m_node[j];

				acc[3*m  ] += (float) t.x;
				acc[3*m+1] += (float) t.y;
				acc[3*m+2] += (float) t.z;
			}

			for (j=0; j<ms.Nodes(); ++j)
			{
				vec3d t = ms.m_Lm[j];// + ss.m_gap[j]*pbi->m_eps;
				int m = ms.m_node[j];

				acc[3*m  ] += (float) t.x;
				acc[3*m+1] += (float) t.y;
				acc[3*m+2] += (float) t.z;
			}
		}

		FESurfaceConstraint* psc = dynamic_cast<FESurfaceConstraint*>(pci);
		if (psc)
		{
			FESurfaceConstraintSurface& ss = psc->m_ss;
			FESurfaceConstraintSurface& ms = psc->m_ms;
			for (j=0; j<ss.Nodes(); ++j)
			{
				vec3d t = ss.m_Lm[j];// + ss.m_gap[j]*pbi->m_eps;
				int m = ss.m_node[j];

				acc[3*m  ] += (float) t.x;
				acc[3*m+1] += (float) t.y;
				acc[3*m+2] += (float) t.z;
			}

			for (j=0; j<ms.Nodes(); ++j)
			{
				vec3d t = ms.m_Lm[j];// + ss.m_gap[j]*pbi->m_eps;
				int m = ms.m_node[j];

				acc[3*m  ] += (float) t.x;
				acc[3*m+1] += (float) t.y;
				acc[3*m+2] += (float) t.z;
			}
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(pci);
		if (pf)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double ti[4], tn[4], gi[4], gn[4], li[4], ln[4];
			int ni, ne;

			int npass = (pf->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FEFacetSlidingSurface& s = (n==0?pf->m_ss:pf->m_ms);

				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k)
					{
						FEFacetSlidingSurface::Data& pt = s.m_Data[j][k];
						li[k] = pt.m_Lm;
						gi[k] = pt.m_gap;
						ti[k] = li[k] + pf->m_epsn*gi[k];

						gi[k] = (gi[k]>=0?gi[k] : 0);
						ti[k] = (ti[k]>=0?ti[k] : 0);
					}

					el.project_to_nodes(li, ln);
					el.project_to_nodes(gi, gn);
					el.project_to_nodes(ti, tn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						acc[3*m  ] += (float) (ln[k]>=0?ln[k]:0);
						acc[3*m+1] += (float) (gn[k]>=0?gn[k]:0);
						acc[3*m+2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.GetMesh().Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[3*j  ] /= (float) val[j]; 
				acc[3*j+1] /= (float) val[j]; 
				acc[3*j+2] /= (float) val[j]; 
			}
		}

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (ps2)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double ti[4], tn[4], gi[4], gn[4], li[4], ln[4];
			int ni, ne;

			int npass = (ps2->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FESlidingSurface2& s = (n==0?ps2->m_ss:ps2->m_ms);

				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k)
					{
						FESlidingSurface2::Data& pt = s.m_Data[j][k];
						li[k] = pt.m_Lmd;
						gi[k] = pt.m_gap;
						ti[k] = li[k] + ps2->m_epsn*gi[k];

						gi[k] = (gi[k]>=0?gi[k] : 0);
						ti[k] = (ti[k]>=0?ti[k] : 0);
					}

					el.project_to_nodes(li, ln);
					el.project_to_nodes(gi, gn);
					el.project_to_nodes(ti, tn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						acc[3*m  ] += (float) (ln[k]>=0?ln[k]:0);
						acc[3*m+1] += (float) (gn[k]>=0?gn[k]:0);
						acc[3*m+2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.GetMesh().Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[3*j  ] /= (float) val[j]; 
				acc[3*j+1] /= (float) val[j]; 
				acc[3*j+2] /= (float) val[j]; 
			}
		}

		FESlidingInterface3* ps3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (ps3)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double ti[4], tn[4], gi[4], gn[4], li[4], ln[4];
			int ni, ne;
			
			int npass = (ps3->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FESlidingSurface3& s = (n==0?ps3->m_ss:ps3->m_ms);
				
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k)
					{
						FESlidingSurface3::Data& pt = s.m_Data[j][k];

						li[k] = pt.m_Lmd;
						gi[k] = pt.m_gap;
						ti[k] = li[k] + ps3->m_epsn*gi[k];
						
						gi[k] = (gi[k]>=0?gi[k] : 0);
						ti[k] = (ti[k]>=0?ti[k] : 0);
					}
					
					el.project_to_nodes(li, ln);
					el.project_to_nodes(gi, gn);
					el.project_to_nodes(ti, tn);
					
					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						acc[3*m  ] += (float) (ln[k]>=0?ln[k]:0);
						acc[3*m+1] += (float) (gn[k]>=0?gn[k]:0);
						acc[3*m+2] += (float) (tn[k]>=0?tn[k]:0);
						val[m]++;
					}
				}
			}
			
			for (j=0; j<fem.GetMesh().Nodes(); ++j) if (val[j] > 1) 
			{ 
				acc[3*j  ] /= (float) val[j]; 
				acc[3*j+1] /= (float) val[j]; 
				acc[3*j+2] /= (float) val[j]; 
			}
		}
	}

	fwrite(&acc[0], sizeof(float)*3, fem.GetMesh().Nodes(), m_fp);
}

//-----------------------------------------------------------------------------
void LSDYNAPlotFile::write_fluid_pressures()
{
	FEMesh& mesh = m_pfem->GetMesh();
	float t;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		t = (float) node.m_pt;
		fwrite(&t, sizeof(float),1, m_fp);
	}
}

//-----------------------------------------------------------------------------
void LSDYNAPlotFile::write_temperatures()
{
	FEMesh& mesh = m_pfem->GetMesh();
	float t;
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);
		t = (float) node.m_T;
		fwrite(&t, sizeof(float), 1, m_fp);
	}
}

//-----------------------------------------------------------------------------
void LSDYNAPlotFile::write_heat_flux()
{
	FEMesh& mesh = m_pfem->GetMesh();

	int i, j;

	vector<vec3d> qn(mesh.Nodes());
	vector<int> val(mesh.Nodes());
	for (i=0; i<mesh.Nodes(); ++i) val[i] = 0;

	vec3d ew;
	int n;
	for (int nd=0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);

				// calculate average flux
				ew = vec3d(0,0,0);
				for (j=0; j<el.GaussPoints(); ++j)
				{
					FEMaterialPoint& mp = *el.m_State[j];
					FEHeatMaterialPoint* pt = (mp.ExtractData<FEHeatMaterialPoint>());

					ew += pt->m_q;
				}

				ew /= el.GaussPoints();
	
				// project to nodes
				for (j=0; j<el.Nodes(); ++j)
				{
					n = el.m_node[j];
					qn[n] += ew;
					++val[n];
				}
			}
		}
	}
	for (i=0; i<mesh.Nodes(); ++i) if (val[i] != 0) qn[i] /= val[i];

	// output nodal fluxes
	float af[3];
	for (i=0; i<mesh.Nodes(); ++i)
	{
		af[0] = (float) qn[i].x;
		af[1] = (float) qn[i].y;
		af[2] = (float) qn[i].z;

		fwrite(af, sizeof(float), 3, m_fp);
	}
}

//-----------------------------------------------------------------------------

void LSDYNAPlotFile::write_contact_pressures()
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	vector<float> t; t.assign(mesh.Nodes(), 0);

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FEContactInterface* pci = fem.ContactInterface(i);

		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(pci);
		if (psi)
		{
			FESlidingSurface& ms = psi->m_ms;
			FESlidingSurface& ss = psi->m_ss;

			double tk;
			double scale = psi->Penalty();
			for (j=0; j<ms.Nodes(); ++j)
			{
				double eps = ms.m_eps[j]*scale;
				tk = ms.m_Lm[j] + eps*ms.m_gap[j];
				tk = MBRACKET(tk);
				t[ms.m_node[j]] += (float) tk;
			}

			for (j=0; j<ss.Nodes(); ++j)
			{
				double eps = ss.m_eps[j]*scale;
				tk = ss.m_Lm[j] + eps*ss.m_gap[j];
				tk = MBRACKET(tk);
				t[ss.m_node[j]] += (float) tk;
			}
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(pci);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.m_node[j]] += (float) ms.m_Lm[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) ss.m_Lm[j].norm();
		}

		FEPeriodicBoundary* pbi = dynamic_cast<FEPeriodicBoundary*>(pci);
		if (pbi)
		{
			FEPeriodicSurface& ms = pbi->m_ms;
			FEPeriodicSurface& ss = pbi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.m_node[j]] += (float) ms.m_Lm[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) ss.m_Lm[j].norm();
		}

		FESurfaceConstraint* psc = dynamic_cast<FESurfaceConstraint*>(pci);
		if (psc)
		{
			FESurfaceConstraintSurface& ms = psc->m_ms;
			FESurfaceConstraintSurface& ss = psc->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.m_node[j]] += (float) ms.m_Lm[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) ss.m_Lm[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(pci);
		if (pri)
		{
			FERigidWallSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) ss.Lm[j];
		}

		FEFacet2FacetSliding* pfs = dynamic_cast<FEFacet2FacetSliding*>(pci);
		if (pfs)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double Li[4], Ln[4];
			int ni, ne, n, k;

			int npass = (pfs->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FEFacetSlidingSurface& s = (n==0?pfs->m_ss:pfs->m_ms);

				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k) 
					{
						FEFacetSlidingSurface::Data& pt = s.m_Data[j][k];
						double eps = pfs->m_epsn*pt.m_eps;
						Li[k] = pt.m_Lm + eps*pt.m_gap;
						Li[k] =  MBRACKET(Li[k]);
					}

					el.project_to_nodes(Li, Ln);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						t[m] += (float) MBRACKET(Ln[k]);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.GetMesh().Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}
	}

	fwrite(&t[0], sizeof(float), fem.GetMesh().Nodes(), m_fp);
}

void LSDYNAPlotFile::write_contact_gaps()
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();

	vector<float> t; t.assign(mesh.Nodes(), 0);

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FEContactInterface* pci = fem.ContactInterface(i);

		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(pci);
		if (psi)
		{
			FESlidingSurface& ms = psi->m_ms;
			FESlidingSurface& ss = psi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.m_node[j]] += (float) (ms.m_gap[j] < 0 ? 0 : ms.m_gap[j]);
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) (ss.m_gap[j] < 0 ? 0 : ss.m_gap[j]);
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(pci);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.m_node[j]] += (float) ms.m_gap[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) ss.m_gap[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(pci);
		if (pri)
		{
			FERigidWallSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.m_node[j]] += (float) (ss.gap[j] < 0? 0 : ss.gap[j]);
		}

		FEFacet2FacetSliding* pf = dynamic_cast<FEFacet2FacetSliding*>(pci);
		if (pf)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double gi[4], gn[4];
			int ni, ne, n, k;

			int npass = (pf->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FEFacetSlidingSurface& s = (n==0?pf->m_ss:pf->m_ms);

				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k)
					{
						double gk = s.m_Data[j][k].m_gap;
						gi[k] = (gk < 0 ? 0 : gk);
					}

					el.project_to_nodes(gi, gn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						double gk = gn[k];
						t[m] += (float) (gk < 0 ? 0 : gk);
						val[m]++;
					}
				}
			}

			for (j=0; j<fem.GetMesh().Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}

		FESlidingInterface2* ps2 = dynamic_cast<FESlidingInterface2*>(pci);
		if (ps2)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double gi[4], gn[4];
			int ni, ne, n, k;

			int npass = (ps2->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FESlidingSurface2& s = (n==0?ps2->m_ss:ps2->m_ms);

				int nint = 0;
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k, ++nint)
					{
						gi[k] = s.m_Data[j][k].m_gap;
					}

					el.project_to_nodes(gi, gn);

					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						t[m] += (float) gn[k];
						val[m]++;
					}
				}
			}

			for (j=0; j<mesh.Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}

		FESlidingInterface3* ps3 = dynamic_cast<FESlidingInterface3*>(pci);
		if (ps3)
		{
			vector<int> val; val.assign(fem.GetMesh().Nodes(), 0);
			double gi[4], gn[4];
			int ni, ne, n, k;
			
			int npass = (ps3->m_btwo_pass?2:1);
			for (n=0; n<npass; ++n)
			{
				FESlidingSurface3& s = (n==0?ps3->m_ss:ps3->m_ms);
				
				for (j=0; j<s.Elements(); ++j)
				{
					FESurfaceElement& el = s.Element(j);
					ne = el.Nodes();
					ni = el.GaussPoints();
					for (k=0; k<ni; ++k)
					{
						gi[k] = s.m_Data[j][k].m_gap;
					}
					
					el.project_to_nodes(gi, gn);
					
					for (k=0; k<ne; ++k)
					{
						int m = el.m_node[k];
						t[m] += (float) gn[k];
						val[m]++;
					}
				}
			}
			
			for (j=0; j<mesh.Nodes(); ++j) if (val[j] > 1) t[j] /= (float) val[j];
		}
	}

	fwrite(&t[0], sizeof(float), fem.GetMesh().Nodes(), m_fp);
}

void LSDYNAPlotFile::write_reaction_forces()
{
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();
	FEAnalysisStep* pstep = dynamic_cast<FEAnalysisStep*>(fem.GetCurrentStep());
	FESolidSolver& solver = dynamic_cast<FESolidSolver&>(*pstep->m_psolver);
	vector<double>& Fr = solver.m_Fr;

	int N = mesh.Nodes(), i;
	vector<float> R(3*N);

	for (i=0; i<N; ++i)
	{
		FENode& node = mesh.Node(i);
		int* id = node.m_ID;
		R[3*i  ] = (float) (-id[0] - 2 >= 0 ? Fr[-id[0]-2] : 0);
		R[3*i+1] = (float) (-id[1] - 2 >= 0 ? Fr[-id[1]-2] : 0);
		R[3*i+2] = (float) (-id[2] - 2 >= 0 ? Fr[-id[2]-2] : 0);
	}

	fwrite(&R[0], sizeof(float)*3, N, m_fp);
}

void LSDYNAPlotFile::write_material_fibers()
{
	int i, j, n, nd;
	FEModel& fem = *m_pfem;
	FEMesh& mesh = fem.GetMesh();
	int N = mesh.Nodes();
	vector<vec3d> v(N);
	for (i=0; i<N; ++i) v[i] = vec3d(0,0,0);

	vec3d r;
	for (nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain* pbd = dynamic_cast<FESolidDomain*>(&mesh.Domain(nd));
		if (pbd)
		{
			for (i=0; i<pbd->Elements(); ++i)
			{
				FESolidElement& el = pbd->Element(i);
				n = el.GaussPoints();
				r = vec3d(0,0,0);
				for (j=0; j<n; ++j)
				{
					FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
					r.x += pt.Q[0][0];
					r.y += pt.Q[1][0];
					r.z += pt.Q[2][0];
				}
				r /= n;

				n = el.Nodes();
				for (j=0; j<n; ++j) v[el.m_node[j]] += r;
			}
		}

		FEShellDomain* psd = dynamic_cast<FEShellDomain*>(&mesh.Domain(nd));
		if (psd)
		{
			for (i=0; i<psd->Elements(); ++i)
			{
				FEShellElement& el = psd->Element(i);
				n = el.GaussPoints();
				r = vec3d(0,0,0);
				for (j=0; j<n; ++j)
				{
					FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();
					r.x += pt.Q[0][0];
					r.y += pt.Q[1][0];
					r.z += pt.Q[2][0];
				}
				r /= n;

				n = el.Nodes();
				for (j=0; j<n; ++j) v[el.m_node[j]] += r;
			}
		}
	}

	float vf[3];
	for (i=0; i<N; ++i)
	{
		v[i].unit();
		vf[0] = (float) v[i].x;
		vf[1] = (float) v[i].y;
		vf[2] = (float) v[i].z;

		fwrite(vf, sizeof(float)*3, 1, m_fp);
	}
}

float LSDYNAPlotFile::fiber_strain(FESolidElement &el, int j)
{
	// see if this element belongs to a tranversely isotropic material
	FEMaterial* pm = m_pfem->GetMaterial(el.GetMatID());
	if (dynamic_cast<FETransverselyIsotropic*>(pm->GetElasticMaterial()))
	{
		FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();

		// get the initial fiber direction
		vec3d a0;
		a0.x = pt.Q[0][0];
		a0.y = pt.Q[1][0];
		a0.z = pt.Q[2][0];

		// calculate the current fiber direction
		vec3d a = pt.F*a0;

		return (float) a.norm();
	}
	else return 0.f;
}

float LSDYNAPlotFile::dev_fiber_strain(FESolidElement &el, int j)
{
	// see if this element belongs to a tranversely isotropic material
	FEMaterial* pm = m_pfem->GetMaterial(el.GetMatID());
	if (dynamic_cast<FETransverselyIsotropic*>(pm->GetElasticMaterial()))
	{
		FEElasticMaterialPoint& pt = *el.m_State[j]->ExtractData<FEElasticMaterialPoint>();

		// get the jacobian
		double J = pt.J;
		double Jm13 = pow(J, -1.0/3.0);

		// get the initial fiber direction
		vec3d a0;
		a0.x = pt.Q[0][0];
		a0.y = pt.Q[1][0];
		a0.z = pt.Q[2][0];

		// calculate the current fiber direction
		vec3d a = pt.F*a0;
		return (float) (Jm13*a.unit());
	}
	else return 0.f;
}

//-----------------------------------------------------------------------------
void LSDYNAPlotFile::write_displacements()
{
	FEMesh& mesh = m_pfem->GetMesh();

	float xf[3];
	for (int i=0; i<mesh.Nodes(); ++i)
	{
		FENode& node = mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_rt.x;
		xf[1] = (float) node.m_rt.y;
		xf[2] = (float) node.m_rt.z;

		fwrite(xf, sizeof(float), 3, m_fp);
	}
}
