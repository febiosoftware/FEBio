// restart module
#include "stdafx.h"
#include "fem.h"
#include "FERestartImport.h"

//-----------------------------------------------------------------------------
//!  This routine reads a binary archive that stores a restart point and prepares
//!  the FEM data to be restarted from this point
//!	\param[in] szfile name of the file

bool FEM::Restart(const char* szfile)
{
	// Open the restart file
	FERestartImport file;
	if (file.Load(*this, szfile) == false)
	{
		char szerr[256];
		file.GetErrorMessage(szerr);
		fprintf(stderr, "%s", szerr);
		return false;
	}

	// Open the log file for appending
	if (m_log.append(m_szlog) == false)
	{
		printf("WARNING: Could not reopen log file. A new log file is created\n");
		m_log.open(m_szlog);
		return false;
	}

	// Open the plot file for appending
	if (m_szplot)
	{
		if (m_plot.Append(*this, m_szplot) == false)
		{
			printf("FATAL ERROR: Failed reopening plot database %s\n", m_szplot);
			return false;
		}
	}

	// inform the user from where the problem is restarted
	m_log.printbox(" - R E S T A R T -", "Restarting from time %lg.\n", m_ftime);

	return true;
}

//-----------------------------------------------------------------------------
//!  Reads or writes the current state to/from a binary file
//!  This is used to restart the solution from a saved position
//!  or to create a restart point.
//!  A version number is written to file to make sure the same
//!  format is used for reading and writing.
//! \param[in] ar the archive to which the data is serialized
//! \sa Archive

bool FEM::Serialize(Archive& ar)
{
	int i, j, k, n, m;

	if (ar.IsSaving())
	{
		// --- version number ---
		ar << RSTRTVERSION;

		// --- analysis data ---
		ar << m_Step.size();
		for (i=0; i<m_Step.size(); ++i) m_Step[i].Serialize(ar);

		ar << m_ftime;

		// --- Load Curve Data ---
		ar << LoadCurves();
		for (i=0; i<LoadCurves(); ++i)
		{
			FELoadCurve& lc = *GetLoadCurve(i);
			n = lc.Points();
			ar << n;
			for (j=0; j<n; ++j)
			{
				LOADPOINT& p = lc.LoadPoint(j);
				ar << p.time << p.value;
			}
		}

		// --- Material Data ---
		ar << Materials();
		for (i=0; i<Materials(); ++i)
		{
			FEMaterial* pmat = GetMaterial(i);

			// store the type string
			ar << pmat->GetTypeString();

			// store the name
			ar << pmat->GetName();

			// store all parameters
			auto_ptr<FEParameterList> pl(pmat->GetParameterList());
			int n = pl->Parameters();
			ar << n;
			list<FEParam>::iterator it = pl->first();
			for (int j=0; j<n; ++j, ++it)
			{
				// store the value
				switch (it->m_itype)
				{
				case FE_PARAM_INT    : ar << it->value<int   >(); break;
				case FE_PARAM_BOOL   : ar << it->value<bool  >(); break;
				case FE_PARAM_DOUBLE : ar << it->value<double>(); break;
				case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<double>()[k]; } break;
				case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar << it->pvalue<int   >()[k]; } break;
				default:
					assert(false);
				}
			}

			// not all parameters can be serialized through the parameter lists
			// so we have to save those parameters the hard way

			if (dynamic_cast<FETransverselyIsotropic*>(pmat))
			{
				FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*>(pmat);
				ar << pm->ca0;
				ar << pm->beta;
				ar << pm->l0;
				ar << pm->refl;
			}

			if (dynamic_cast<FERigid*>(pmat))
			{
				FERigid* pm = dynamic_cast<FERigid*>(pmat);
				ar.write(pm->m_bc, sizeof(int), 6);
				ar.write(pm->m_fc, sizeof(int), 6);
				ar.write(pm->m_fs, sizeof(double), 6);
			}
		}

		// --- Geometry Data ---
		m_mesh.Serialize(ar);

		// write solid element state data
		for (i=0; i<m_mesh.SolidElements(); ++i)
		{
			FESolidElement& el = m_mesh.SolidElement(i);
			for (j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}

		// write shell element state data
		for (i=0; i<m_mesh.ShellElements(); ++i)
		{
			FEShellElement& el = m_mesh.ShellElement(i);
			for (j=0; j<el.GaussPoints(); ++j) el.m_State[j]->Serialize(ar);
		}

		// surface elements
		n = m_psurf->Elements();
		ar << n;
		for (i=0; i<n; ++i)
		{
			FESurfaceElement& el = m_psurf->Element(i);
			ar << el.Type();
			ar << el.GetMatID();
			ar << el.m_nID;
			ar << el.m_nrigid;
			ar << el.m_node;
			ar << el.m_lnode;
		}

		// rigid bodies
		ar << m_nreq << m_nrb << m_nrm;
		for (i=0; i<m_nrb; ++i)
		{
			FERigidBody& rb = m_RB[i];

			ar << rb.m_nID << rb.m_mat << rb.m_mass << rb.m_Fr << rb.m_Mr;
			ar << rb.m_r0 << rb.m_rt << rb.m_qt;
			ar.write(rb.m_bc, sizeof(int), 6);
			ar.write(rb.m_LM, sizeof(int), 6);
			ar.write(rb.m_Up, sizeof(double), 6);
			ar.write(rb.m_Ut, sizeof(double), 6);
			ar.write(rb.m_du, sizeof(double), 6);
		}

		// rigid joints
		ar << m_nrj;
		for (i=0; i<m_nrj; ++i)
		{
			ar.write(&m_RJ[i], sizeof(FERigidJoint), 1);
		}

		// contact data
		ar << m_bcontact;
		ar << ContactInterfaces();
		for (i=0; i<ContactInterfaces(); ++i)
		{
			FESlidingInterface* psi = dynamic_cast<FESlidingInterface*>(&m_CI[i]);

			if (psi)
			{
				ar << psi->Type();
				ar << psi->npass;
				ar << psi->m_eps;
				ar << psi->m_atol;
				ar << psi->m_nplc;
				ar << psi->nse;
				ar << psi->nme;
				ar << psi->m_bautopen;
				ar << psi->m_nsegup;
				ar << psi->m_blaugon;

				for (j=0; j<2; ++j)
				{
					FEContactSurface& s = (j==0? psi->m_ss : psi->m_ms);

					int ne = s.Elements();
					ar << ne;

					for (k=0; k<ne; ++k)
					{
						FESurfaceElement& el = s.Element(k);
						ar << el.Type();
						ar << el.GetMatID() << el.m_nID << el.m_nrigid;
						ar << el.m_node;
						ar << el.m_lnode;
					}

					ar << s.eps;
					ar << s.gap;
					ar << s.nu;
					ar << s.rs;
					ar << s.Lm;
					ar << s.off;
				}
			}

			FETiedInterface* pti = dynamic_cast<FETiedInterface*>(&m_CI[i]);
			if (pti)
			{
				ar << pti->Type();
				ar << pti->m_eps;
				ar << pti->m_atol;
				ar << pti->m_nplc;
				ar << pti->nse;
				ar << pti->nme;

				for (j=0; j<2; ++j)
				{
					FETiedContactSurface& s = (j==0? pti->ss : pti->ms);

					int ne = s.Elements();
					ar << ne;

					for (k=0; k<ne; ++k)
					{
						FESurfaceElement& el = s.Element(k);
						ar << el.Type();
						ar << el.GetMatID() << el.m_nID << el.m_nrigid;
						ar << el.m_node;
						ar << el.m_lnode;
					}

					ar << s.gap;
					ar << s.rs;
					ar << s.Lm;
				}
			}

			FERigidWallInterface* pri = dynamic_cast<FERigidWallInterface*>(&m_CI[i]);
			if (pri)
			{
				ar << pri->Type();
				ar << pri->m_eps;
				ar << pri->m_atol;
				ar << pri->m_nplc;

				FEContactSurface& s = pri->m_ss;

				int ne = s.Elements();
				ar << ne;

				for (k=0; k<ne; ++k)
				{
					FESurfaceElement& el = s.Element(k);
					ar << el.Type();
					ar << el.GetMatID() << el.m_nID << el.m_nrigid;
					ar << el.m_node;
					ar << el.m_lnode;
				}

				ar << s.gap;
				ar << s.nu;
				ar << s.rs;
				ar << s.Lm;
				ar << s.off;
				
				// plane data
				if (dynamic_cast<FEPlane*>(pri->m_mp))
				{
					FEPlane* pp = dynamic_cast<FEPlane*>(pri->m_mp);
					ar << FE_RIGID_PLANE;
					ar << pp->m_nplc;
					double* a = pp->GetEquation();
					ar << a[0] << a[1] << a[2] << a[3];
				}
				else if (dynamic_cast<FERigidSphere*>(pri->m_mp))
				{
					FERigidSphere* ps = dynamic_cast<FERigidSphere*>(pri->m_mp);
					ar << FE_RIGID_SPHERE;
					ar << ps->m_rc;
					ar << ps->m_R;
					ar << ps->m_nplc[0] << ps->m_nplc[1] << ps->m_nplc[2];
				}
			}
		}

		// --- Boundary Condition Data ---

		// displacements
		ar << m_DC.size();
		for (i=0; i<m_DC.size(); ++i)
		{
			ar.write(m_DC+i, sizeof(FENodalDisplacement), 1);
		}

		// nodal loads
		ar << m_FC.size();
		for (i=0; i<m_FC.size(); ++i)
		{
			ar.write(m_FC+i, sizeof(FENodalForce), 1);
		}

		// pressure forces
		ar << m_PC.size();
		for (i=0; i<m_PC.size(); ++i)
		{
			ar.write(m_PC+i, sizeof(FE_FACE_PRESSURE), 1);
		}

		// body forces
		ar.write(m_BF  ,sizeof(FE_BODYFORCE), 1);
		ar.write(m_BF+1,sizeof(FE_BODYFORCE), 1);
		ar.write(m_BF+2,sizeof(FE_BODYFORCE), 1);

		ar << m_acc;

		// discrete elements
		ar << m_DE.size();
		for (i=0; i<m_DE.size(); ++i)
		{
			FE_DISCRETE_ELEMENT& de = m_DE[i];
			ar << de.n1 << de.n2;
			ar << de.E;
		}

		// --- Direct Solver Data ---
		ar << m_nsolver;
		ar << m_neq;

		// --- I/O-stuff ---
		ar << m_szfile << m_szplot << m_szlog << m_szdump;
		ar << m_sztitle;
	}
	else
	{
		int mat=0;

		// --- version ---
		ar >> n;

		// --- analysis data ---
		int nsteps;
		ar >> nsteps;
		for (i=0; i<nsteps; ++i)
		{
			FEAnalysis* pstep = new FEAnalysis(*this);
			m_Step.add(pstep);
			pstep->Serialize(ar);
		}

		ar >> m_ftime;

		// --- Load Curve Data ---
		int nlc;
		ar >> nlc;
		for (i=0; i<nlc; ++i)
		{
			FELoadCurve* plc = new FELoadCurve;
			ar >> n;
			plc->Create(n);
			for (j=0; j<n; ++j)
			{
				LOADPOINT& p = plc->LoadPoint(j);
				ar >> p.time >> p.value;
			}
			AddLoadCurve(plc);
		}

		// --- Material Data ---
		int nmat;
		char szmat[256] = {0}, szvar[256] = {0};
		ar >> nmat;
		for (i=0; i<nmat; ++i)
		{
			// read the type string
			ar >> szmat;

			// create a material
			FEMaterial* pmat = FEMaterialFactory::CreateMaterial(szmat);
			assert(pmat);
			AddMaterial(pmat);

			// read the name
			ar >> szmat;
			pmat->SetName(szmat);

			// read all parameters
			auto_ptr<FEParameterList> pl(pmat->GetParameterList());
			int n = 0;
			ar >> n;
			assert(n == pl->Parameters());
			list<FEParam>::iterator it = pl->first();
			for (int j=0; j<n; ++j, ++it)
			{
				// read the value
				switch (it->m_itype)
				{
				case FE_PARAM_INT    : ar >> it->value<int   >(); break;
				case FE_PARAM_BOOL   : ar >> it->value<bool  >(); break;
				case FE_PARAM_DOUBLE : ar >> it->value<double>(); break;
				case FE_PARAM_DOUBLEV: { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<double>()[k]; } break;
				case FE_PARAM_INTV   : { for (int k=0; k<it->m_ndim; ++k) ar >> it->pvalue<int   >()[k]; } break;
				default:
					assert(false);
				}
			}
	
			// not all parameters can be serialized through the parameter lists
			// so we have to save those parameters the hard way

			if (dynamic_cast<FETransverselyIsotropic*>(pmat))
			{
				FETransverselyIsotropic* pm = dynamic_cast<FETransverselyIsotropic*>(pmat);
				ar >> pm->ca0;
				ar >> pm->beta;
				ar >> pm->l0;
				ar >> pm->refl;
			}

			if (dynamic_cast<FERigid*>(pmat))
			{
				FERigid* pm = dynamic_cast<FERigid*>(pmat);
				ar.read(pm->m_bc, sizeof(int), 6);
				ar.read(pm->m_fc, sizeof(int), 6);
				ar.read(pm->m_fs, sizeof(double), 6);
			}
		}

		// --- Geometry Data ---

		m_mesh.Serialize(ar);

		// read solid element state data
		for (i=0; i<m_mesh.SolidElements(); ++i)
		{
			FESolidElement& el = m_mesh.SolidElement(i);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}

		// read shell element state data
		for (i=0; i<m_mesh.ShellElements(); ++i)
		{
			FEShellElement& el = m_mesh.ShellElement(i);
			for (j=0; j<el.GaussPoints(); ++j)
			{
				el.SetMaterialPointData(GetMaterial(el.GetMatID())->CreateMaterialPointData(), j);
				el.m_State[j]->Serialize(ar);
			}
		}

		// surface elements
		ar >> n;
		if (n) 
		{
			m_psurf->Create(n);

			for (i=0; i<n; ++i)
			{
				FESurfaceElement& el = m_psurf->Element(i);
				ar >> m;
				el.SetType(m);

				ar >> mat; el.SetMatID(mat);
				ar >> el.m_nID;
				ar >> el.m_nrigid;
				ar >> el.m_node;
				ar >> el.m_lnode;
			}

			// initialize surface data
			m_psurf->Init();
		}

		// rigid bodies
		ar >> m_nreq >> m_nrb >> m_nrm;
		if (m_nrb) m_RB.create(m_nrb);
		for (i=0; i<m_nrb; ++i)
		{
			FERigidBody& rb = m_RB[i];

			ar >> rb.m_nID >> rb.m_mat >> rb.m_mass >> rb.m_Fr >> rb.m_Mr;
			ar >> rb.m_r0 >> rb.m_rt >> rb.m_qt;
			ar.read(rb.m_bc, sizeof(int), 6);
			ar.read(rb.m_LM, sizeof(int), 6);
			ar.read(rb.m_Up, sizeof(double), 6);
			ar.read(rb.m_Ut, sizeof(double), 6);
			ar.read(rb.m_du, sizeof(double), 6);

			rb.AttachToFEM(this);
		}

		// rigid joints
		ar >> m_nrj;
		for (i=0; i<m_nrj; ++i)
		{
			FERigidJoint* prj = new FERigidJoint(this);
			ar.read(prj, sizeof(FERigidJoint), 1);
			m_RJ.add(prj);
		}

		// contact data
		int numci;
		ar >> m_bcontact;
		ar >> numci;
		for (i=0; i<numci; ++i)
		{
			ar >> n;

			switch (n)
			{
			case FE_CONTACT_SLIDING:
				{
					FESlidingInterface* ps = new FESlidingInterface(this);
					FESlidingInterface& si = *ps;
					m_CI.add(ps);

					ar >> si.npass;
					ar >> si.m_eps;
					ar >> si.m_atol;
					ar >> si.m_nplc;
					ar >> si.nse;
					ar >> si.nme;
					ar >> si.m_bautopen;
					ar >> si.m_nsegup;
					ar >> si.m_blaugon;

					if (si.m_nplc >= 0) si.m_pplc = &m_LC[si.m_nplc];

					for (j=0; j<2; ++j)
					{
						FEContactSurface& s = (j==0? si.m_ss : si.m_ms);

						int ne=0;
						ar >> ne;
						s.Create(ne);

						for (k=0; k<ne; ++k)
						{
							FESurfaceElement& el = s.Element(k);
	
							ar >> n;
							el.SetType(n);
		
							ar >> mat >> el.m_nID >> el.m_nrigid;
							ar >> el.m_node;
							ar >> el.m_lnode;

							el.SetMatID(mat);
						}

						// initialize surface
						s.Init();

						// read the contact data
						// Note that we do this after Init() since this data gets 
						// initialized to zero there
						ar >> s.eps;
						ar >> s.gap;
						ar >> s.nu;
						ar >> s.rs;
						ar >> s.Lm;
						ar >> s.off;
					}
				}
				break;

			case FE_CONTACT_TIED:
				{
					FETiedInterface* pt = new FETiedInterface(this);
					FETiedInterface& ti = *pt;
					m_CI.add(pt);

					ar >> ti.m_eps;
					ar >> ti.m_atol;
					ar >> ti.m_nplc;
					ar >> ti.nse;
					ar >> ti.nme;

					if (ti.m_nplc >= 0) ti.m_pplc = &m_LC[ti.m_nplc];

					for (j=0; j<2; ++j)
					{
						FETiedContactSurface& s = (j==0? ti.ss : ti.ms);

						int ne=0;
						ar >> ne;
						s.Create(ne);

						for (k=0; k<ne; ++k)
						{
							FESurfaceElement& el = s.Element(k);
	
							ar >> n;
							el.SetType(n);
		
							ar >> mat >> el.m_nID >> el.m_nrigid;
							ar >> el.m_node;
							ar >> el.m_lnode;

							el.SetMatID(mat);
						}

						// initialize surface
						s.Init();

						// read the contact data
						// Note that we do this after Init() since this data gets 
						// initialized to zero there
						ar >> s.gap;
						ar >> s.rs;
						ar >> s.Lm;
					}
				}
				break;

			case FE_CONTACT_RIGIDWALL:
				{
					FERigidWallInterface* pr = new FERigidWallInterface(this);
					FERigidWallInterface& ri = *pr;
					m_CI.add(pr);

					ar >> ri.m_eps;
					ar >> ri.m_atol;
					ar >> ri.m_nplc;


					if (ri.m_nplc >= 0) ri.m_pplc = &m_LC[ri.m_nplc];

					FEContactSurface& s = ri.m_ss;

					int ne=0;
					ar >> ne;
					s.Create(ne);

					for (k=0; k<ne; ++k)
					{
						FESurfaceElement& el = s.Element(k);
	
						ar >> n;
						el.SetType(n);
		
						ar >> mat >> el.m_nID >> el.m_nrigid;
						ar >> el.m_node;
						ar >> el.m_lnode;

						el.SetMatID(mat);
					}

					// initialize surface
					s.Init();

					ar >> s.gap;
					ar >> s.nu;
					ar >> s.rs;
					ar >> s.Lm;
					ar >> s.off;

					// plane data
					int ntype;
					ar >> ntype;
					switch (ntype)
					{
					case FE_RIGID_PLANE:
						{
							ri.SetMasterSurface(new FEPlane(this));
							FEPlane& pl = dynamic_cast<FEPlane&>(*ri.m_mp);
							ar >> pl.m_nplc;
							if (pl.m_nplc >= 0) pl.m_pplc = &m_LC[pl.m_nplc];
							double* a = pl.GetEquation();
							ar >> a[0] >> a[1] >> a[2] >> a[4];
						}
						break;
					case FE_RIGID_SPHERE:
						{
							ri.SetMasterSurface(new FERigidSphere(this));
							FERigidSphere& s = dynamic_cast<FERigidSphere&>(*ri.m_mp);
							ar >> s.m_rc;
							ar >> s.m_R;
							ar >> s.m_nplc[0] >> s.m_nplc[1] >> s.m_nplc[2];
						}
						break;
					default:
						assert(false);
					}
				}
			}
		}

		// --- Boundary Condition Data ---

		// displacements
		int ndis;
		ar >> ndis;
		m_DC.create(ndis);
		for (i=0; i<ndis; ++i)
		{
			ar.read(m_DC+i, sizeof(FENodalDisplacement), 1);
		}

		// nodal loads
		int ncnf;
		ar >> ncnf;
		m_FC.create(ncnf);
		for (i=0; i<ncnf; ++i)
		{
			ar.read(m_FC+i, sizeof(FENodalForce), 1);
		}

		// pressure forces
		int npr;
		ar >> npr;
		m_PC.create(npr);
		for (i=0; i<npr; ++i)
		{
			ar.read(m_PC+i, sizeof(FE_FACE_PRESSURE), 1);
		}

		// body forces
		ar.read(m_BF  ,sizeof(FE_BODYFORCE), 1);
		ar.read(m_BF+1,sizeof(FE_BODYFORCE), 1);
		ar.read(m_BF+2,sizeof(FE_BODYFORCE), 1);

		ar >> m_acc;

		// discrete elements
		int nde;
		ar >> nde;
		if (nde > 0)
		{
			m_DE.setsize(nde);
			for (i=0; i<nde; ++i)
			{
				FE_DISCRETE_ELEMENT& de = m_DE[i];
				ar >> de.n1 >> de.n2;
				ar >> de.E;
			}
		}

		// --- Direct Solver Data ---
		ar >> m_nsolver;
		ar >> m_neq;

		// --- I/O-stuff ---
		ar >> m_szfile >> m_szplot >> m_szlog >> m_szdump;
		ar >> m_sztitle;
	}

	return true;
}
