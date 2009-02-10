// PlotFile.cpp: implementation of the PlotFile class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "PlotFile.h"
#include "fem.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlotFile::PlotFile()
{
	m_bsstrn = false;
	m_nfield[0] = -1;	// displacement
	m_nfield[1] = -1;	// velocity
	m_nfield[2] = -1;	// acceleration
	m_nfield[3] = -1;	// temperature
}

PlotFile::~PlotFile()
{
	// close the archive
	m_ar.Close();
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : PlotFile::Open
//  Opens the PLOT database. That is, an archive is created and the PLOTHEADER 
//  is written as well as the initial geometry of the model.
//

bool PlotFile::Open(FEM& fem, const char* szfile)
{
	int i, j, N;

	// open the archive
	if (m_ar.Create(szfile) == false) return false;

	m_pfem = &fem;

	int itype = fem.m_pStep->m_itype;

	// check the field values
	if (m_nfield[0] == -1) m_nfield[0] = PLOT_DISPLACEMENT;
	if (m_nfield[1] == -1)
	{
		m_nfield[1] = PLOT_NONE;
		if ((itype == FE_STATIC_PORO) || (itype == FE_DYNAMIC)) m_nfield[1] = PLOT_VELOCITY;
	}
	if (m_nfield[2] == -1)
	{
		m_nfield[2] = PLOT_NONE;
		if (itype == FE_STATIC_PORO) m_nfield[2] = PLOT_FLUID_FLUX;
		else if (itype == FE_DYNAMIC) m_nfield[2] = PLOT_ACCELERATION;
		else if (fem.m_bcontact) m_nfield[2] = PLOT_CONTACT_TRACTION;
	}
	if (m_nfield[3] == -1)
	{
		m_nfield[3] = PLOT_NONE;
		if (itype == FE_STATIC_PORO) m_nfield[3] = PLOT_FLUID_PRESSURE;
		else if (fem.m_bcontact) m_nfield[3] = PLOT_CONTACT_GAP;
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
	plh.nel2  = 0;
	plh.nel4  = fem.m_mesh.ShellElements();
	plh.nel8  = fem.m_mesh.SolidElements();
	plh.nglbv = 0;
	plh.nummat2 = 0;
	plh.nummat4 = 0;
	plh.nummat8 = fem.Materials();
	plh.nump    = fem.m_mesh.Nodes();
	plh.nv1d    = 6;
	plh.nv2d    = (m_bsstrn? 44 : 32);
	plh.nv3d    = 7;
	plh.maxint  = 0;
	plh.neips   = 2000;
	plh.ioshl4  = 0;

	// store the header
	m_ph = plh;

	m_ar.write(&plh, sizeof(PLOTHEADER), 1);

	// write the material coordinates
	float xf[3];
	for (i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		xf[0] = (float) node.m_r0.x;
		xf[1] = (float) node.m_r0.y;
		xf[2] = (float) node.m_r0.z;

		m_ar.write(xf, sizeof(float), 3);
	}

	// write the connectivity and material number
	// Note that we increment all numbers by 1 since
	// the plot database expects 1-based arrays
	int n[9];

	FEMesh& mesh = fem.m_mesh;

	int nid = 1;

	// write solid element data
	// note that we reindex all elements so that the ID
	// corresponds to the nr in the plot file
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);

		el.m_nID = nid++;
		
		N = el.Nodes();
		switch (el.Type())
		{
		case FE_HEX:
		case FE_RIHEX:
			for (j=0; j<N; ++j) n[j] = el.m_node[j]+1;
			break;
		case FE_PENTA:
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
		case FE_TET:
			n[0] = el.m_node[0]+1;
			n[1] = el.m_node[1]+1;
			n[2] = el.m_node[2]+1;
			n[3] = el.m_node[2]+1;
			n[4] = n[5] = n[6] = n[7] = el.m_node[3]+1;
			break;
		}
			
		n[8] = el.GetMatID()+1;
		
		m_ar.write(n, sizeof(int), 9);
	}

	// write shell element data
	for (i=0; i<mesh.ShellElements(); ++i)
	{
		FEShellElement& el = mesh.ShellElement(i);

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
		
		m_ar.write(n, sizeof(int), 5);
	}

	return true;
}


///////////////////////////////////////////////////////////////////////////////
// FUNCTION : PlotFile::Append
// Opens a PLOT database for appending
//
// TODO: I always need to make sure that the same amount of data is written
// to the file after opening for appending. In other words, make sure that
// the flags are the same as when opening the file for the first time.

bool PlotFile::Append(FEM& fem, const char* szfile)
{
	// write the header
	PLOTHEADER plh = {0};

	// open the file
	m_ar.Open(szfile);

	// read the plot file header
	m_ar.read(&plh, sizeof(PLOTHEADER), 1);

	if (plh.nv2d == 32) m_bsstrn = false;
	else if (plh.nv2d == 44) m_bsstrn = true;
	else 
	{
		m_ar.Close();
		return false;
	}

	// close the file
	m_ar.Close();

	m_pfem = &fem;

	// reopen the plot file for appending
	return (m_ar.Append(szfile));
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : PlotFile::Write
//  Writes the current time step data to the PLOT database. The data archive
//  must have been opened with the Open or Append function.
//

bool PlotFile::Write(FEM& fem)
{
	// make sure the archive is opened
	if (!m_ar.IsValid()) return false;

	// Ok, let's proceed
	int i, j;

	// write the time value
	float time = (float) fem.m_ftime;
	m_ar << time;

	FEMesh& mesh = fem.m_mesh;

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
		case PLOT_CONTACT_TRACTION: write_contact_tractions();
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
		default:
			assert(false);
		}
	}

	// write the stresses
	// note that at this point the plastic strain field (= s[6]) is not written to
	// or better, is always set to zero.

	// write solid element data
	float s[44] = {0};
	double f;
	int nint, neln;
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);

		for (j=0; j<7; ++j) s[j] = 0;

		nint = el.GaussPoints();

		f = 1.0 / (double) nint;

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		// we output the average stress values of the gauss points
		for (j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint& pt = *(el.m_State[j]->ExtractData<FEElasticMaterialPoint>());

			s[0] += (float) (f*pt.s.xx());
			s[1] += (float) (f*pt.s.yy());
			s[2] += (float) (f*pt.s.zz());
			s[3] += (float) (f*pt.s.xy());
			s[4] += (float) (f*pt.s.yz());
			s[5] += (float) (f*pt.s.xz());
		}

		m_ar.write(s, sizeof(float), 7);
	}

	// write shell element data
	mat3ds E;
	for (i=0; i<mesh.ShellElements(); ++i)
	{
		for (j=0; j<44; ++j) s[j] = 0;

		FEShellElement& el = mesh.ShellElement(i);
		if (!el.IsRigid())
		{
			try
			{
				mesh.UnpackElement(el);

				neln = el.Nodes();
				f = 1.0 / (double) neln;

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
		m_ar.write(s, sizeof(float), (m_bsstrn?44:32));
	}

	// now flush the archive to make sure we don't loose any results
	m_ar.Flush();

	return true;
}

//-----------------------------------------------------------------------------

void PlotFile::Close()
{
	m_ar.Close();
}

//-----------------------------------------------------------------------------

void PlotFile::write_displacements()
{
	FEM& fem = *m_pfem;

	float xf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		xf[0] = (float) node.m_rt.x;
		xf[1] = (float) node.m_rt.y;
		xf[2] = (float) node.m_rt.z;

		m_ar.write(xf, sizeof(float), 3);
	}
}

void PlotFile::write_velocities()
{
	FEM& fem = *m_pfem;
	float vf[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		vf[0] = (float) node.m_vt.x;
		vf[1] = (float) node.m_vt.y;
		vf[2] = (float) node.m_vt.z;

		m_ar.write(vf, sizeof(float), 3);
	}
}

void PlotFile::write_accelerations()
{
	FEM& fem = *m_pfem;
	float af[3];
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);

		// since the PLOT file requires floats we need to convert
		// the doubles to single precision
		af[0] = (float) node.m_at.x;
		af[1] = (float) node.m_at.y;
		af[2] = (float) node.m_at.z;

		m_ar.write(af, sizeof(float), 3);
	}
}

void PlotFile::write_fluid_flux()
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	int i, j;

	vector<vec3d> wn(fem.m_mesh.Nodes());
	vector<int> val(fem.m_mesh.Nodes());
	for (i=0; i<fem.m_mesh.Nodes(); ++i) val[i] = 0;

	vec3d ew;
	int n;
	for (i=0; i<mesh.SolidElements(); ++i)
	{
		FESolidElement& el = mesh.SolidElement(i);

		// calculate average flux
		ew = vec3d(0,0,0);
		for (j=0; j<el.GaussPoints(); ++j) 
		{
			FEMaterialPoint& mp = *el.m_State[j];
			FEPoroElasticMaterialPoint* pt = (mp.ExtractData<FEPoroElasticMaterialPoint>());

			ew += pt->m_w;
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
	for (i=0; i<fem.m_mesh.Nodes(); ++i) if (val[i] != 0) wn[i] /= val[i];

	// output nodal fluxes
	float af[3];
	for (i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		af[0] = (float) wn[i].x;
		af[1] = (float) wn[i].y;
		af[2] = (float) wn[i].z;

		m_ar.write(af, sizeof(float), 3);
	}
}

void PlotFile::write_contact_tractions()
{
	FEM& fem = *m_pfem;

	int i, j;

	vector<float[3]> acc(fem.m_mesh.Nodes());
	for (i=0; i<fem.m_mesh.Nodes(); ++i) acc[i][0] = acc[i][1] = acc[i][2] = 0;
	for (i=0; i<fem.m_CI.size(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*> (&fem.m_CI[i]);
		if (psi)
		{
			FEContactSurface& ss = psi->m_ss;
			for (j=0; j<ss.Nodes(); ++j)
			{
				int m = ss.node[j];
				acc[m][0] = (float) ss.Lt[j][0];
				acc[m][1] = (float) ss.Lt[j][1];
				acc[m][2] = (float) ss.Lm[j];
			}
		}
	}

	m_ar.write(acc, sizeof(float)*3, fem.m_mesh.Nodes() );
}

void PlotFile::write_fluid_pressures()
{
	FEM& fem = *m_pfem;
	float t;
	for (int i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		FENode& node = fem.m_mesh.Node(i);
		t = (float) node.m_pt;
		m_ar << t;
	}
}

void PlotFile::write_contact_pressures()
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	vector<float> t(mesh.Nodes());
	t.zero();

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(&fem.m_CI[i]);
		if (psi)
		{
			FEContactSurface& ms = psi->m_ms;
			FEContactSurface& ss = psi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.Lm[j];
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.Lm[j];
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(&fem.m_CI[i]);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.Lm[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.Lm[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&fem.m_CI[i]);
		if (pri)
		{
			FEContactSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.Lm[j];
		}
	}
	
	for (i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		m_ar << t[i];
	}
}

void PlotFile::write_contact_gaps()
{
	FEM& fem = *m_pfem;
	FEMesh& mesh = fem.m_mesh;

	vector<float> t(mesh.Nodes());
	t.zero();

	int i, j;

	for (i=0; i<fem.ContactInterfaces(); ++i)
	{
		FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(&fem.m_CI[i]);
		if (psi)
		{
			FEContactSurface& ms = psi->m_ms;
			FEContactSurface& ss = psi->m_ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) (ms.gap[j] < 0 ? 0 : ms.gap[j]);
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) (ss.gap[j] < 0 ? 0 : ss.gap[j]);
		}

		FETiedInterface* pti = dynamic_cast<FETiedInterface*>(&fem.m_CI[i]);
		if (pti)
		{
			FETiedContactSurface& ms = pti->ms;
			FETiedContactSurface& ss = pti->ss;

			for (j=0; j<ms.Nodes(); ++j) t[ms.node[j]] += (float) ms.gap[j].norm();
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) ss.gap[j].norm();
		}

		FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&fem.m_CI[i]);
		if (pri)
		{
			FEContactSurface& ss = pri->m_ss;
			for (j=0; j<ss.Nodes(); ++j) t[ss.node[j]] += (float) (ss.gap[j] < 0? 0 : ss.gap[j]);
		}
	}
	
	for (i=0; i<fem.m_mesh.Nodes(); ++i)
	{
		m_ar << t[i];
	}
}
