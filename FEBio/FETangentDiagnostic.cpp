// FETangentDiagnostic.cpp: implementation of the FETangentDiagnostic class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "FETangentDiagnostic.h"
#include "FEBox.h"

void print_matrix(Logfile& log, matrix& m)
{
	int i, j;
	int N = m.rows();
	int M = m.columns();

	log.printf("\n    ");
	for (i=0; i<N; ++i) log.printf("%15d ", i);
	log.printf("\n----");
	for (i=0; i<N; ++i) log.printf("----------------", i);

	for (i=0; i<N; ++i)
	{
		log.printf("\n%2d: ", i);
		for (j=0; j<M; ++j)
		{
			log.printf("%15lg ", m[i][j]);
		}
	}
	log.printf("\n");
}

void mnbrak(double* pa, double* pb, double* pc, double* pfa, double* pfb, double* pfc, double (*func)(double));
double golden(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FETangentDiagnostic* FETangentDiagnostic::m_pthis = 0;


FETangentDiagnostic::FETangentDiagnostic(FEM& fem) : FEDiagnostic(fem)
{
	m_strain = 0;
}

FETangentDiagnostic::~FETangentDiagnostic()
{

}

bool FETangentDiagnostic::Init()
{
	FEMesh& mesh = m_fem.m_mesh;

	// assign the material to the mesh
	mesh.SetMatID(0);

	// TODO: find an easier way to create material point data
	for (int i=0; i<8; ++i)
	{
		mesh.SolidElement(0).SetMaterialPointData( m_fem.GetMaterial(0)->CreateMaterialPointData(), i);
	}

	return FEDiagnostic::Init();
}

bool FETangentDiagnostic::Run()
{
	FEM& fem = m_fem;

	FEMesh& mesh = fem.m_mesh;

	Logfile& log = fem.m_log;

	Logfile::MODE oldmode = log.SetMode(Logfile::FILE_ONLY);

	// get and initialize the first step
	fem.m_Step[0].Init();

	// get and initialize the solver
	FESolver& solver = *fem.m_pStep->m_psolver;
	solver.Init();

	// solve the problem
	solve();

	// get the one and only element
	FESolidElement& el = mesh.SolidElement(0);
	mesh.UnpackElement(el);

	// set up the element stiffness matrix
	matrix k0(24, 24);
	k0.zero();
	solver.ElementStiffness(el, k0);

	// print the element stiffness matrix
	log.printf("\nActual stiffness matrix:\n");
	print_matrix(log, k0);

	// now calculate the derivative of the residual
	matrix k1;
	deriv_residual(k1);

	// print the approximate element stiffness matrix
	log.printf("\nApproximate stiffness matrix:\n");
	print_matrix(log, k1);

	// finally calculate the difference matrix
	log.printf("\n");
	matrix kd(24, 24);
	double kmax = 0, kij;
	int i0 = -1, j0 = -1, i, j;
	for (i=0; i<24; ++i)
		for (j=0; j<24; ++j)
		{
			kd[i][j] = k0[i][j] - k1[i][j];
			kij = 100.0*fabs(kd[i][j] / k0[i][j]);
			if (kij > kmax) 
			{
				kmax = kij;
				i0 = i;
				j0 = j;
			}
		}

	// print the difference
	log.printf("\ndifference matrix:\n");
	print_matrix(log, kd);

	log.SetMode(oldmode);

	log.printf("\nMaximum difference: %lg%% (at (%d,%d))\n", kmax, i0, j0);

	return (kmax < 1e-3);
}

void FETangentDiagnostic::deriv_residual(matrix& ke)
{
	FEM& fem = m_fem;

	// get the solver
	FESolver& solver = *fem.m_pStep->m_psolver;

	// get the mesh
	FEMesh& mesh = fem.m_mesh;

	// get the one and only element
	FESolidElement& el = mesh.SolidElement(0);
	mesh.UnpackElement(el);

	// first calculate the initial residual
	vector<double> f0(24);
	f0.zero();
	solver.InternalForces(el, f0);

	// now calculate the perturbed residuals
	ke.Create(24, 24);
	ke.zero();
	int i, j, nj;
	int N = mesh.Nodes();
	double dx = 1e-8;
	vector<double> f1(24);
	for (j=0; j<3*N; ++j)
	{
		FENode& node = mesh.Node(el.m_node[j/3]);
		nj = j%3;

		switch (nj)
		{
		case 0: node.m_rt.x += dx; break;
		case 1: node.m_rt.y += dx; break;
		case 2: node.m_rt.z += dx; break;
		}


		solver.UpdateStresses();
		mesh.UnpackElement(el);

		f1.zero();
		solver.InternalForces(el, f1);

		switch (nj)
		{
		case 0: node.m_rt.x -= dx; break;
		case 1: node.m_rt.y -= dx; break;
		case 2: node.m_rt.z -= dx; break;
		}

		solver.UpdateStresses();
		mesh.UnpackElement(el);

		for (i=0; i<3*N; ++i) ke[i][j] = -(f1[i] - f0[i])/dx;
	}
}

double FETangentDiagnostic::residual(double d)
{
	FETangentDiagnostic& dia = *m_pthis;

	FEM& fem = dia.m_fem;
	FEMesh& mesh = fem.m_mesh;
	FESolver& solver = *fem.m_pStep->m_psolver;

	FESolidElement& el = mesh.SolidElement(0);
	mesh.UnpackElement(el);

	// set the deformation
	for (int n=0; n<8; ++n)
	{
		FENode& node = mesh.Node(n);
		node.m_rt.x = node.m_r0.x*(1+dia.m_stretch);
		node.m_rt.y = node.m_r0.y*(1 + d);
		node.m_rt.z = node.m_r0.z*(1 + d);
	}

	// update stresses
	solver.UpdateStresses();

	// get the residual
	vector<double> R(24);
	R.zero();
	solver.InternalForces(el, R);

	double r = sqrt(R*R);

	return r;
}

void FETangentDiagnostic::solve()
{
	FEM& fem = m_fem;
	FEMesh& mesh = fem.m_mesh;

	// calculate stretch
	m_stretch = sqrt(2.0*m_strain + 1) - 1;

	// set the static this pointer
	m_pthis = this;

	double a = -0.046;
	double b = -0.045;
	double c, fa, fb, fc;
	mnbrak(&a, &b, &c, &fa, &fb, &fc, FETangentDiagnostic::residual);

	double beta = 0;
	double res = golden(a, b, c, FETangentDiagnostic::residual, 1e-5, &beta);

	double d = residual(beta);

	if (fem.m_pStep->m_nplot != FE_PLOT_NEVER) fem.m_plot.Write(fem);
}

bool FETangentDiagnostic::ParseSection(XMLTag& tag)
{
	if (tag != "Scenario") return false;

	FEM& fem = m_fem;

	const char* szname = tag.AttributeValue("name");
	if (strcmp(szname, "uni-axial") == 0)
	{
		// create the mesh
		FEBox box;
		box.Create(1, 1, 1, vec3d(0,0,0), vec3d(1,1,1));

		// set the geometry
		fem.m_mesh = box;

		++tag;
		while (!tag.isend())
		{
			if (tag == "strain") tag.value(m_strain);
			else throw XMLReader::InvalidTag(tag);
			++tag;
		}
	}
	else throw XMLReader::InvalidAttributeValue(tag, "name", szname);

	return true;
}

#define SHFT2(a, b, c) (a)=(b);(b)=(c);
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a, b) ((b)>=0?fabs(a):(-fabs(a)))
#define FMAX(a, b) ((a)>(b)?(a):(b))

void mnbrak(double* pa, double* pb, double* pc, double* pfa, double* pfb, double* pfc, double (*func)(double))
{
	const double GOLD = 1.618034;
	const double TINY = 1.0e-20;
	const double GLIMIT = 100;

	double& a = *pa;
	double& b = *pb;
	double& c = *pc;

	double& fa = *pfa;
	double& fb = *pfb;
	double& fc = *pfc;

	double ulim, u, r, q, fu, dum;

	fa = func(a);
	fb = func(b);
	if (fb>fa)
	{
		SHFT(dum, a, b, dum);
		SHFT(dum, fb, fa, dum);
	}

	c = b+GOLD*(b - a);
	fc = func(c);
	while (fb > fc)
	{
		r = (b - a)*(fb - fc);
		q = (b - c)*(fb - fa);
		u = b - ((b - c)*q - (b - a)*r) / (2.0*SIGN(FMAX(fabs(q-r), TINY), q-r));

		ulim = b + GLIMIT*(c - b);

		if ((b - u)*(u - c) > 0)
		{
			fu = func(u);
			if (fu < fc)
			{
				a = b;
				b = u;
				fa = fb;
				fb = fu;
				return;
			}
			else if (fu > fb)
			{
				c = u;
				fc = fu;
				return;
			}

			u = c + GOLD*(c - b);
			fu = func(u);
		}
		else if ((c - u)*(u - ulim) > 0)
		{
			fu = func(u);
			if (fu < fc)
			{
				SHFT(b, c, u, c + GOLD*(c - b));
				SHFT(fb, fc, fu, func(u));
			}
		}
		else if ((u-ulim)*(ulim - c) >= 0)
		{
			u = ulim;
			fu = func(u);
		}
		else
		{
			u = c + GOLD*(c - b);
			fu = func(u);
		}

		SHFT(a, b, c, u);
		SHFT(fa, fb, fc, fu);
	}
}

double golden(double ax, double bx, double cx, double (*f)(double), double tol, double* xmin)
{
	const double R = 0.61803399;
	const double C = 1 - R;

	double f1, f2, x0, x1, x2, x3;

	x0 = ax;
	x3 = cx;

	if (fabs(cx - bx) > fabs(bx - ax))
	{
		x1 = bx;
		x2 = bx + C*(cx - bx);
	}
	else
	{
		x2 = bx;
		x1 = bx - C*(bx - ax);
	}
	f1 = f(x1);
	f2 = f(x2);

	while (fabs(x3 - x0) > tol*(fabs(x1) + fabs(x2)))
	{
		if (f2 < f1)
		{
			SHFT(x0, x1, x2, R*x1 + C*x3);
			SHFT2(f1, f2, f(x2));
		}
		else
		{
			SHFT(x3, x2, x1, R*x2 + C*x0);
			SHFT2(f2, f1, f(x1));
		}
	}

	if (f1 < f2)
	{
		*xmin = x1;
		return f1;
	}
	else
	{
		*xmin = x2;
		return f2;
	}
}
