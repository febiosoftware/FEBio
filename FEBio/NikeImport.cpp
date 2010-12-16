// NikeImport.cpp: implementation of the FENIKEImport class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "NikeImport.h"
#include "FENeoHookean.h"
#include "FEMooneyRivlin.h"
#include "FETransIsoMooneyRivlin.h"
#include "FERigid.h"
#include "FESolidSolver.h"
#include "FESlidingInterface.h"

#include <string.h>
#include <stdlib.h>

#define ABS(x) ((x)>=0?(x):(-(x)))

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : read_line
// This function reads a line from a text file, while skipping over comment
// lines (ie. lines which begin with an asterisk (*))
// The function returns 0 if an end-of-file or an error is encountered.
// Otherwise it returns szline;
//

char* FENIKEImport::read_line(FILE* fp, char* szline, int n)
{
	// read a line while skipping over comment lines
	while (fgets(szline, n, fp) && (szline[0]=='*'));

	// make sure we did not encounter a problem
	if (feof(fp) || ferror(fp)) return NULL;

	// remove newline
	char* ch = strrchr(szline, '\n');
	if (ch) *ch = 0;

	return szline;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::Load
// Imports a NIKE input deck
//

bool FENIKEImport::Load(FEM& fem, const char* szfile)
{
	// open the file
	if (Open(szfile, "rb") == 0) return errf("FATAL ERROR: Failed opening input file %s\n\n", szfile);

	// keep a copy of the input file name
	fem.SetInputFilename(szfile);

	// make sure we have a solver defined
	if (fem.m_pStep->m_psolver == 0) fem.m_pStep->m_psolver = new FESolidSolver(fem);

	// Read the control deck
	if (ReadControlDeck(fem) == false) return false;

	// Read the material deck
	if (ReadMaterialDeck(fem) == false) return false;

	// Read the geometry
	if (ReadGeometry(fem) == false) return false;

	// Read the load curves
	if (ReadCurveDeck(fem) == false) return false;

	// Read boundary conditions and body force decks
	if (ReadBCDecks(fem) == false) return false;

	// close the input file
	Close();

	// all done!
	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::ReadControlDeck
//
// Reads the loadcurve deck from a nike input file
//

bool FENIKEImport::ReadControlDeck(FEM& fem)
{
	// error constants
	const int ERR_EOF		= 0;
	const int ERR_CC		= 1;
	const int ERR_FORMAT	= 2;
	const int ERR_ATYPE		= 3;

	// error message formats
	const char szerr[][256] = {
		"FATAL ERROR: Unexpected end of file in input file %s\n\n",
		"FATAL ERROR: Error reading control card %d\n\n",
		"FATAL ERROR: Incorrect input format. Must be FL\n\n",
		"FATAL ERROR: Unrecognized analysis type on CC 7\n\n"
	};

	// storage for input line
	const int MAX_LINE = 256;
	char szline[MAX_LINE];
	int nread;

	// -------- control card 1 --------
	// read title
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

	// -------- control card 2 --------
	// read number of materials, nodes and elements
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
	char szf[5] = {0};
	nread = sscanf(szline, "%2s%3d%10d%10d%*10d%10d%*5d%5d%*5d%*5d%5d", szf, &m_nmat, &m_nn, &m_nbel, &m_nsel, &m_numsi, &m_nrn);
	if (nread != 7) return errf(szerr[ERR_CC], 2);

	// make sure input format is correct
	if (strcmp(szf,"FL") != 0) return errf(szerr[ERR_FORMAT]);

	// -------- control card 3 --------
	// read nr of timesteps, time step size, auto step flag, opt nr, of iterations, min step, max step, ...
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
	char szauto[5] = {0};
	nread = sscanf(szline, "%10d%10lg %4s%5d%5d%10lg%10lg", &fem.m_pStep->m_ntime, &fem.m_pStep->m_dt0, szauto, &fem.m_pStep->m_maxretries, &fem.m_pStep->m_iteopt, &fem.m_pStep->m_dtmin, &fem.m_pStep->m_dtmax);
	if (nread != 7) return errf(szerr[ERR_CC], 3);
	if ((strcmp(szauto,"auto") == 0) || (strcmp(szauto,"AUTO") == 0)) fem.m_pStep->m_bautostep = true;

	if (fem.m_pStep->m_iteopt == 0) fem.m_pStep->m_iteopt = 11;
	if (fem.m_pStep->m_maxretries == 0) fem.m_pStep->m_maxretries = 5;
	if (fem.m_pStep->m_dtmin == 0) fem.m_pStep->m_dtmin = fem.m_pStep->m_dt0 / 3.0;

	if (fem.m_pStep->m_dtmax == 0) fem.m_pStep->m_dtmax = fem.m_pStep->m_dt0 * 3.0;
	else if (fem.m_pStep->m_dtmax < 0) fem.m_pStep->m_nmplc = (int) (-fem.m_pStep->m_dtmax) - 1;

	// -------- control card 4 --------
	// read nr of loadcurves, nr of displ bound cards, nr of press bound cards, nr of conc. nodal force cards
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
	int ax, ay, az;
	nread = sscanf(szline, "%5d%*5d%5d%5d%5d%*5d%*5d%5d%5d%5d", &m_nlc, &m_ncnf, &m_npr, &m_ndis, &ax, &ay, &az);
	if (nread != 7) return errf(szerr[ERR_CC], 4);

	if ((ax!=0)||(ay!=0)||(az!=0))
	{
		FEBodyForce* pbf = new FEBodyForce;
		if (ax != 0) pbf->lc[0] = 1; else pbf->lc[0] = -1;
		if (ay != 0) pbf->lc[1] = 1; else pbf->lc[1] = -1;
		if (az != 0) pbf->lc[2] = 1; else pbf->lc[2] = -1;
		fem.m_BF.push_back(pbf);
	}

	// -------- control card 5 --------
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
	int irfreq;
	nread = sscanf(szline, "%*5d%*5d%*5d%*5d%*5d%*5d%*5d%5d", &irfreq);
	if (nread != 1) return errf(szerr[ERR_CC], 5);

	fem.m_pStep->m_bDump = (irfreq != 0);

	// -------- control card 6 --------
	// read convergence tolerances and max nr of iteration values
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
	double Dtol, Etol, Rtol, LStol;
	int maxups, maxref;
	int sops;
	nread = sscanf(szline,"%*5d%5d%*10d%*10d%5d%5d%10lg%10lg%10lg%10lg", &sops, &maxups, &maxref, &Dtol, &Etol, &Rtol, &LStol);
	if (nread != 7) return errf(szerr[ERR_CC], 6);

	if (sops != 0) fem.m_bwopt = 1;

	if (maxups ==   0) maxups = 10;
	if (maxref ==   0) maxref = 15;
	if (Dtol   == 0.0) Dtol   = 0.001;
	if (Etol   == 0.0) Etol   = 0.01;
	if (Rtol   == 0.0) Rtol   = 0; // i.e. deactivated
	if (LStol  == 0.0) LStol  = 0.9;

	fem.m_pStep->m_psolver->m_Dtol = Dtol;
	fem.m_pStep->m_psolver->m_Etol = Etol;
	fem.m_pStep->m_psolver->m_Rtol = Rtol;
	fem.m_pStep->m_psolver->m_LStol = LStol;
	fem.m_pStep->m_psolver->m_bfgs.m_maxups = maxups;
	fem.m_pStep->m_psolver->m_bfgs.m_maxref = maxref;

	// -------- control card 7 --------
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

	// -------- control card 8 --------
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
//	nread = sscanf(szline, "%5d", &fem.m_istiffpr);

	// -------- control card 9 --------
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

	// -------- control card 10 --------
	if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::ReadMaterialDeck
//
// Reads the material deck from a nike input file
//

bool FENIKEImport::ReadMaterialDeck(FEM& fem)
{
	int i, j, m;

	// error constants
	const int ERR_EOF		= 0;
	const int ERR_MAT		= 1;
	const int ERR_NOMAT		= 2;
	const int ERR_AOPT      = 3;

	// error message formats
	const char szerr[][256] = {
		"FATAL ERROR: Unexpected end of file in input file %s\n\n",
		"FATAL ERROR: Error in material card %d of material %d\n\n",
		"FATAL ERROR: Material %d is not supported\n\n",
		"FATAL ERROR: Invalid value for material axis option on CC6 of material %d\n\n"
	};

	const int MAX_LINE = 256;
	char szline[MAX_LINE];
	char szname[MAX_LINE];
	int nread;

	int mat, itype;
	double mp[8][8] = {0}, density;
	for (i=0; i<m_nmat; ++i)
	{
		// -------- material card 1 --------
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		nread = sscanf(szline, "%*5d%5d%10lg%5d", &mat, &density, &itype);
		if (nread != 3) return errf(szerr[ERR_MAT], 1, i);

		m = 6;
		if (itype == 2) m = 8;

		// -------- material card 2 --------
		if (read_line(m_fp, szname, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		// -------- material card 3 - 8 --------
		for (j=0; j<m; j++)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

			nread = sscanf(szline, "%lg%lg%lg%lg%lg%lg%lg%lg", mp[j],mp[j]+1,mp[j]+2,mp[j]+3,mp[j]+4,mp[j]+5,mp[j]+6,mp[j]+7);
			if (nread != 8) return errf(szerr[ERR_MAT], j+3, i+1);
		}

		// check material type
		switch (mat)
		{
		case 1:	// Currently mapped to a neo-Hookean
			{
				double E, v;
				FENeoHookean* pmat = new FENeoHookean;
				pmat->SetName(szname);
				fem.AddMaterial(pmat);
				E = mp[0][0];
				v = mp[1][0];

				pmat->m_E = E;
				pmat->m_v = v;
				pmat->m_density = density;
			}
			break;

		case 15: // Mooney-Rivlin
			{
				FEMooneyRivlin* pmat = new FEMooneyRivlin;
				pmat->SetName(szname);
				fem.AddMaterial(pmat);

				double A, B, v;
				A = mp[0][0];
				B = mp[1][0];
				v = mp[2][0];

				pmat->c1 = A;
				pmat->c2 = B;

				pmat->m_K = 4.0*(A+B)*(1.0+v)/(3.0*(1.0 - 2.0*v));

				pmat->m_density = density;
			}
			break;

		case 18: // Transversely Isotropic Mooney-Rivlin
			{
				FETransIsoMooneyRivlin* pmat = new FETransIsoMooneyRivlin;
				pmat->SetName(szname);
				fem.AddMaterial(pmat);

				int n1, n2, n4;

				pmat->c1 = mp[0][0];
				pmat->c2 = mp[0][1];
				pmat->m_fib.m_c3 = mp[0][2];
				pmat->m_fib.m_c4 = mp[0][3];
				pmat->m_fib.m_c5 = mp[0][4];

				pmat->m_K  = mp[1][0];
				pmat->m_fib.m_lam1 = mp[1][1];

				int naopt = (int) mp[3][0];

				switch (naopt)
				{
				case FE_FIBER_LOCAL:
					{
						FELocalMap* pmap = new FELocalMap();
						pmat->m_pmap = pmap;

						n1 = (int) mp[4][0];
						n2 = (int) mp[4][1];
						n4 = (int) mp[4][2];

						if ((n1==0) && (n2==0) && (n4==0)) { n1 = 1; n2 = 2; n4 = 4; }

						pmap->SetLocalNodes(n1-1, n2-1, n4-1);
					}
					break;
				case FE_FIBER_SPHERICAL:
					{
						FESphericalMap* pmap = new FESphericalMap();
						pmat->m_pmap = pmap;

						vec3d c;

						c.x = mp[4][0];
						c.y = mp[4][1];
						c.z = mp[4][2];

						pmap->SetSphereCenter(c);
					}
					break;
				case FE_FIBER_VECTOR:
					{
						FEVectorMap* pmap = new FEVectorMap();
						pmat->m_pmap = pmap;

						vec3d a, d;

						a.x = mp[4][0];
						a.y = mp[4][1];
					    a.z = mp[4][2];

						d.x = mp[5][0];
						d.y = mp[5][1];
					    d.z = mp[5][2];

						pmap->SetVectors(a, d);
					}
					break;
				default:
					return errf(szerr[ERR_AOPT], i+1);
				}

				pmat->m_density = density;

				// active contraction input
				pmat->m_fib.m_lcna = (int) mp[5][3] - 1;
				pmat->m_fib.m_ca0  = mp[5][4];
				pmat->m_fib.m_beta = mp[5][5];
				pmat->m_fib.m_refl = mp[5][6];
				pmat->m_fib.m_l0   = mp[5][7];
			}
			break;

		case 20: // Rigid body
			{
				double E, v;
				FERigidMaterial* pmat = new FERigidMaterial;
				pmat->SetName(szname);
				fem.AddMaterial(pmat);
				E = mp[0][0];
				v = mp[1][0];

				pmat->m_E = E;
				pmat->m_v = v;
				pmat->m_density = density;

				pmat->m_com = (int) mp[3][0];
				pmat->m_rc.x = mp[3][1];
				pmat->m_rc.y = mp[3][2];
				pmat->m_rc.z = mp[3][3];

				pmat->m_bc[0] = (int) mp[2][0];
				pmat->m_bc[1] = (int) mp[2][1];
				pmat->m_bc[2] = (int) mp[2][2];
				pmat->m_bc[3] = (int) mp[2][3];
				pmat->m_bc[4] = (int) mp[2][4];
				pmat->m_bc[5] = (int) mp[2][5];

				pmat->m_fc[0] = -1; pmat->m_fs[0] = 1;
				pmat->m_fc[1] = -1; pmat->m_fs[1] = 1;
				pmat->m_fc[2] = -1; pmat->m_fs[2] = 1;
				pmat->m_fc[3] = -1; pmat->m_fs[3] = 1;
				pmat->m_fc[4] = -1; pmat->m_fs[4] = 1;
				pmat->m_fc[5] = -1; pmat->m_fs[5] = 1;
			}
			break;

		default:
			// Oh, oh, we didn't recognize the material type or don't support it.
			return errf(szerr[ERR_NOMAT], mat);
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::ReadGeometry
//
// Reads the nodal and element decks from a nike input file
//

bool FENIKEImport::ReadGeometry(FEM& fem)
{
	int i, j;

	// error constants
	const int ERR_EOF		= 0;
	const int ERR_NODE		= 1;
	const int ERR_HEX		= 2;
	const int ERR_RB		= 3;
	const int ERR_SI		= 4;
	const int ERR_SIT		= 5;
	const int ERR_SHELL		= 6;

	// error message formats
	const char szerr[][256] = {
		"FATAL ERROR: Unexpected end of file in input file %s\n\n",
		"FATAL ERROR: Error in nodal coordinates deck at node %d\n\n",
		"FATAL ERROR: Error in hexahedral deck at element %d\n\n",
		"FATAL ERROR: Error in rigd node deck. Material %d is not a rigid body\n\n",
		"FATAL ERROR: Error encounted readin sliding surface deck\n\n",
		"FATAL ERROR: Error in sliding surface deck. Only interface types 3 are supported\n\n",
		"FATAL ERROR: Error in shell element deck at element %d\n\n"
	};

	const int MAX_LINE = 256;
	char szline[MAX_LINE];
	int nread;

	// create mesh
	FEMesh& mesh = fem.m_mesh;
	mesh.CreateNodes(m_nn);

	FEElasticSolidDomain* pbd = (m_nbel>0 ? new FEElasticSolidDomain(&mesh, 0) : 0);
	FEElasticShellDomain* psd = (m_nsel>0 ? new FEElasticShellDomain(&mesh, 0) : 0);

	if (pbd) { pbd->create(m_nbel); mesh.AddDomain(pbd); }
	if (psd) { psd->create(m_nsel); mesh.AddDomain(psd); }

	int GID[2], nc=0;
	if (pbd) GID[0] = nc++; else GID[0] = -1;
	if (psd) GID[1] = nc++; else GID[1] = -1;

	////////////////////// N O D A L   C O O R D I N A T E S   D E C K //////////////////////

	// initialize ID array
	for (i=0; i<m_nn; ++i)
		for (j=0; j<MAX_NDOFS; ++j) fem.m_mesh.Node(i).m_ID[j] = -1;

	// read geometry
	int bcu, rc;
	double r[3] = {0};
	for (i=0; i<m_nn; ++i)
	{
		// read initial coordinates
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		rc = 0;
		nread = sscanf(szline, "%*8d%5d%20lg%20lg%20lg%5d", &bcu, r, r+1, r+2, &rc);
		if (nread < 4) return errf(szerr[ERR_NODE], i+1);

		FENode& node = fem.m_mesh.Node(i);
		node.m_rt = node.m_r0 = vec3d(r[0], r[1], r[2]);

		// set nodal displacement constraints
		// bc = 1 : constrained x-displacement
		// bc = 2 : constrained y-displacement
		// bc = 3 : constrained z-displacement
		// bc = 4 : constrained x,y-displacement
		// bc = 5 : constrained y,z-displacement
		// bc = 6 : constrained x,z-displacement
		// bc = 7 : constrained x,y,z-displacement
		if ((bcu==0) || (bcu==2) || (bcu==3) || (bcu==5)) node.m_ID[0] = 0;
		if ((bcu==0) || (bcu==1) || (bcu==3) || (bcu==6)) node.m_ID[1] = 0;
		if ((bcu==0) || (bcu==1) || (bcu==2) || (bcu==4)) node.m_ID[2] = 0;

		// set rotational constraints
		if ((rc==0) || (rc==2) || (rc==3) || (rc==5)) node.m_ID[3] = 0;
		if ((rc==0) || (rc==1) || (rc==3) || (rc==6)) node.m_ID[4] = 0;
		if ((rc==0) || (rc==1) || (rc==2) || (rc==4)) node.m_ID[5] = 0;
	}

	//////////////////////////// H E X   E L E M E N T   D E C K ////////////////////////////
	int n[9];
	for (i=0; i<m_nbel; ++i)
	{
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		nread = sscanf(szline,"%*8d%5d%8d%8d%8d%8d%8d%8d%8d%8d", n+8,n,n+1,n+2,n+3,n+4,n+5,n+6,n+7);
		if (nread != 9) return errf(szerr[ERR_HEX], i+1);

		// see what type of element we are dealing with
		assert(pbd);
		FESolidElement& el = pbd->Element(i);

		el.m_gid = 0;

		// since arrays in C are zero-based we need do decrease the nodes and material number
		if ((n[7]==n[3]) && (n[6]==n[3]) && (n[5]==n[3]) && (n[4]==n[3]))
		{
			el.SetType(FE_TET);
			el.m_node[0] = n[0] - 1;
			el.m_node[1] = n[1] - 1;
			el.m_node[2] = n[2] - 1;
			el.m_node[3] = n[3] - 1;
		}
		else if ((n[7]==n[5]) && (n[6]==n[5]))
		{
			// note the strange mapping. This is because NIKE's wedge elements use a
			// different node numbering.
			el.SetType(FE_PENTA);
			el.m_node[0] = n[0]-1;
			el.m_node[1] = n[4]-1;
			el.m_node[2] = n[1]-1;
			el.m_node[3] = n[3]-1;
			el.m_node[4] = n[6]-1;
			el.m_node[5] = n[2]-1;
		}
		else
		{
			el.SetType(FE_HEX);
			for (j=0; j<8; j++) el.m_node[j] = n[j] - 1;
		}
		el.SetMatID(n[8] - 1);

		// determine whether this element is rigid
		FERigidMaterial* pm = dynamic_cast<FERigidMaterial*>(fem.GetMaterial(el.GetMatID()));
		if (pm) el.m_nrigid = 0;
	}

	//////////////////////////// S H E L L   E L E M E N T   D E C K ////////////////////////////
	float f[5];
	int N;
	for (i=0; i<m_nsel; ++i)
	{
		// read the first shell card
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		nread = sscanf(szline,"%*8d%5d%8d%8d%8d%8d", n+4,n,n+1,n+2,n+3);
		if (nread != 5) return errf(szerr[ERR_SHELL], i+1);

		// read the second shell card
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		nread = sscanf(szline,"%10g%10g%10g%10g", f, f+1, f+2, f+3);
		if (nread != 4) return errf(szerr[ERR_SHELL], i+1);

		// get the element
		assert(psd);
		FEShellElement& el = psd->Element(i);
		if (n[3] == n[2])
			el.SetType(FE_SHELL_TRI);
		else
			el.SetType(FE_SHELL_QUAD);

		el.m_gid = 1;

		N = el.Nodes();
		for (j=0; j<N; ++j) el.m_node[j] = n[j]-1;
		el.SetMatID(n[4]-1);

		for (j=0; j<N; ++j) el.m_h0[j] = f[j];
	}

	// assign material point data
	for (i=0; i<mesh.Domains(); ++i)
	{
		FEDomain& d = mesh.Domain(i);
		d.InitMaterialPointData();
	}

	////////////////////// R I G I D   N O D E   A N D   F A C E T   D E C K ////////////////

	// initialize array
	for (i=0; i<m_nn; ++i)
	{
		FENode& node = fem.m_mesh.Node(i);
		node.m_rid = -1;
	}

	int nrn = 0;

	int rn[4], rb;
	for (i=0; i<m_nrn; ++i)
	{
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);
		nread = sscanf(szline, "%5d%8d%8d%8d%8d", &rb, rn,rn+1,rn+2,rn+3);

		if (dynamic_cast<FERigidMaterial*>(fem.GetMaterial(rb-1)) == 0)
		{
			return errf(szerr[ERR_RB], rb);
		}

		for (j=0; j<nread-1; ++j, ++N, ++nrn)
		{
			FERigidNode* prn = new FERigidNode;
			prn->nid = rn[j]-1;
			prn->rid = rb-1;
			fem.m_RN.push_back(prn);
		}
	}
	// since we may have read less than 4*m_nrn rigid nodes
	// we need to make sure we have the correct size for m_RN
	fem.m_RN.resize(nrn);

	/////////////////////////// S L I D I N G   S U R F A C E   D E C K /////////////////////

	if (m_numsi)
	{
		// read the sliding interface control cards
		int iaux;
		vector<SLIDING_INTERFACE> SI(m_numsi);
		for (i=0; i<m_numsi; ++i)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_SI]);
			nread = sscanf(szline, "%8d%8d%4d%10lg%10lg%*10lg%*10lg%*10lg%*5d%5d", &SI[i].nss, &SI[i].nms, &SI[i].itype, &SI[i].sfac, &SI[i].mus, &iaux);
			if (nread != 6) return errf(szerr[ERR_SI]);

			// make sure interface type == 3
			if (ABS(SI[i].itype) != 3) return errf(szerr[ERR_SIT]);

			// read the auxiliary control card
			if (iaux)
			{
				if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_SI]);
				nread = sscanf(szline, "%5d%10lg%10lg%10lg", &SI[i].iaug, &SI[i].altoln, &SI[i].altolt, &SI[i].tkmult);
				if (nread != 4) return errf(szerr[ERR_SI]);
			}
		}

		// create the sliding interfaces
		for (i=0; i<m_numsi; ++i)
		{
			// create a new sliding interface and add it to the model
			FESlidingInterface* psi = new FESlidingInterface(&fem);
			fem.m_CI.push_back(psi);
			FESlidingInterface& si = *psi;

			// allocate storage for contact surfaces
			si.m_ss.create(SI[i].nss);
			si.m_ms.create(SI[i].nms);

			// set contact parameters
			si.m_eps = SI[i].sfac;
			si.m_atol = SI[i].altoln;
			si.m_mu = SI[i].mus;

			// validate the penalty factor
			if (si.m_eps == 0) si.m_eps = 1;
			if (si.m_eps < 0) si.m_eps *= -1;

			si.m_npass = 2;
			if (SI[i].itype == -3) si.m_npass = 1;

			// set default lag aug tolerance
			if (si.m_atol == 0) si.m_atol = 0.1;
		}

		// read the sliding interfaces
		int en[4], k, N;
		for (i=0; i<m_numsi; ++i)
		{
			FESlidingInterface& si = dynamic_cast<FESlidingInterface&>(*fem.m_CI[i]);

			int nss = si.m_ss.Elements();
			int nms = si.m_ms.Elements();

			// read the slave facets
			for (j=0; j<nss; ++j)
			{
				FESurfaceElement& el = si.m_ss.Element(j);

				if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_SI]);
				nread = sscanf(szline, "%*8d%8d%8d%8d%8d", en, en+1, en+2, en+3);

				if (en[2] == en[3])
					el.SetType(FE_NITRI);
				else
					el.SetType(FE_NIQUAD);

				N = el.Nodes();

				for (k=0; k<N; ++k) el.m_node[k] = en[k]-1;
			}

			// read the master facets
			for (j=0; j<nms; ++j)
			{
				FESurfaceElement& el = si.m_ms.Element(j);

				if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_SI]);
				nread = sscanf(szline, "%*8d%8d%8d%8d%8d", en, en+1, en+2, en+3);

				if (en[2] == en[3])
					el.SetType(FE_NITRI);
				else
					el.SetType(FE_NIQUAD);

				N = el.Nodes();

				for (k=0; k<N; ++k) el.m_node[k] = en[k]-1;
			}
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::ReadCurveDeck
//
// Reads the load curve deck from a nike input file
//

bool FENIKEImport::ReadCurveDeck(FEM& fem)
{
	int i, j;

	// error constants
	const int ERR_EOF		= 0;
	const int ERR_LC		= 1;

	// error message formats
	const char szerr[][256] = {
		"FATAL ERROR: Unexpected end of file in input file %s\n\n",
		"FATAL ERROR: Error in load curve deck of load curve %d\n\n",
	};

	const int MAX_LINE = 256;
	char szline[MAX_LINE];
	int nread;

	int lcs; // size of loadcurve
	double time, val;
	for (i=0; i<m_nlc; ++i)
	{
		// -------- load card 1 --------
		if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

		nread = sscanf(szline, "%*5d%5d", &lcs);
		if (nread != 1) return errf(szerr[ERR_LC], i+1);

		FELoadCurve* plc = new FELoadCurve();
		plc->Create(lcs);

		// -------- load card 2 - n --------
		for (j=0; j<lcs; ++j)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

			nread = sscanf(szline, "%10lg%10lg", &time, &val);
			if (nread != 2) return errf(szerr[ERR_LC], i+1);

			plc->SetPoint(j, time, val);
		}

		fem.AddLoadCurve(plc);
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
// FUNCTION : FENIKEImport::ReadBCDecks
//
// Reads the decks describing condcentrated nodal forces, prescribed displacements,
// pressure forces and body forced from a nike input file
//

bool FENIKEImport::ReadBCDecks(FEM& fem)
{
	int i;

	// error constants
	const int ERR_EOF		= 0;
	const int ERR_CFORCE	= 1;
	const int ERR_DISP		= 2;
	const int ERR_PRESS		= 3;
	const int ERR_BFORCE	= 4;
	const int ERR_DISP_BC	= 5;
	const int ERR_CFORCE_BC	= 6;

	// error message formats
	const char szerr[][256] = {
		"FATAL ERROR: Unexpected end of file in input file %s\n\n",
		"FATAL ERROR: Error in concentrated nodal force deck\n\n",
		"FATAL ERROR: Error in displacement boundary deck\n\n",
		"FATAL ERROR: Error in pressure boundary condition deck\n\n",
		"FATAL ERROR: Error in base acceleration body force deck\n\n",
		"FATAL ERROR: Invalid boundary condition in displacement deck\n\n",
		"FATAL ERROR: Invalid boundary condition in nodal force deck\n\n"
	};

	const int MAX_LINE = 256;
	char szline[MAX_LINE];
	int nread;

	/////////////// C O N C E N T R A T E D   N O D A L   F O R C E   D E C K ///////////////
	if (m_ncnf>0)
	{
		int bc;
		for (i=0; i<m_ncnf; ++i)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

			FENodalForce* pfc = new FENodalForce;

			nread = sscanf(szline, "%8d%5d%5d%10lg", &(pfc->node), &bc, &(pfc->lc), &(pfc->s));
			if (nread != 4) return errf(szerr[ERR_CFORCE]);

			if ((bc < 1) || (bc > 3)) return errf(szerr[ERR_CFORCE_BC]);

			// make indices zero-based
			pfc->node--;	// node point number
			pfc->bc = bc-1;	// direction of load
			pfc->lc--;	// load curve number

			// set default scale factor if requested
			if (pfc->s == 0.0) pfc->s = 1.0;
			fem.m_FC.push_back(pfc);
		}
	}

	////////////////////// P R E S S U R E   B O U N D A R Y   D E C K //////////////////////
	if (m_npr > 0)
	{
		fem.m_psurf = new FEPressureSurface(&fem.m_mesh);
		FEPressureSurface& ps = *fem.m_psurf;
		ps.create(m_npr);
		int n[4], N, j;
		for (i=0; i<m_npr; ++i)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

			// create a surface element
			FESurfaceElement& el = ps.Element(i);

			// read the card data
			FEPressureLoad& pc = ps.PressureLoad(i);
			double* s = pc.s;
			pc.face = i;

			nread = sscanf(szline, "%5d%8d%8d%8d%8d%10lg%10lg%10lg%10lg", &(pc.lc),n,n+1,n+2,n+3,s,s+1,s+2,s+3);
			if (nread != 9) return errf(szerr[ERR_PRESS]);

			if (n[3] == n[2])
				el.SetType(FE_TRI);
			else
				el.SetType(FE_QUAD);

			N = el.Nodes();
			for (j=0; j<N; ++j) el.m_node[j] = n[j]-1;

			pc.lc--;

			// we set the material to -1 to make sure this element is not associated with any material
			el.SetMatID(-1);

			// if all scale factors are zero, default them to 1
			if ((s[0]==0) && (s[1]==0) && (s[2]==0) && (s[3]==0)) s[0] = s[1] = s[2] = s[3] = 1.0;
		}
	}

	////////////////// D I S P L A C E M E N T   B O U N D A R Y   D E C K //////////////////
	if (m_ndis > 0)
	{
		fem.m_DC;
		for (i=0; i<m_ndis; ++i)
		{
			if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

			FENodalDisplacement* pdc = new FENodalDisplacement;
			fem.m_DC.push_back(pdc);
			FENodalDisplacement& dc = *pdc;

			nread = sscanf(szline, "%8d%5d%5d%10lg", &(dc.node), &(dc.bc), &(dc.lc), &(dc.s));
			if (nread != 4) return errf(szerr[ERR_DISP]);

			if ((dc.bc < 1) || (dc.bc > 3)) return errf(szerr[ERR_DISP_BC]);

			// make indices zero-based
			dc.node--;	// node point number
			dc.bc--;	// direction of displacement
			dc.lc--;	// load curve number

			FENode& node = fem.m_mesh.Node(dc.node);

			// make sure that the prescribed dofs are free
			node.m_ID[dc.bc] = 0;
		}
	}

	/////////// B A S E   A C C E L E R A T I O N   B O D Y   F O R C E   D E C K ///////////
	if (!fem.m_BF.empty())
	{
		FEBodyForce& BF = *fem.m_BF[0];
		for (i=0; i<3; ++i)
		{
			// the card is only present if the according flag on CC4 was not zero
			if (BF.lc[i] >= 0)
			{
				if (read_line(m_fp, szline, MAX_LINE) == NULL) return errf(szerr[ERR_EOF], m_szfile);

				nread = sscanf(szline, "%5d%10lg", &(BF.lc[i]), &(BF.s[i]));
				if (nread != 2) return errf(szerr[ERR_BFORCE]);
				BF.lc[i]--;
				if (BF.s[i] == 0) BF.s[i] = 1.0;
			}
		}
	}

	return true;
}
