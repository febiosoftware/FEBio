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
#include <stdio.h>
#include <stdlib.h>
#include "PardisoProjectSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>

//! This implementation of the Pardiso solver is for the version
//! available in the Intel MKL.


#ifdef PARDISODL

#include <thread>

#ifdef PARDISODL_RUNTIME
//-----------------------------------------------------------------------------
// Run-time binding of the Panua Pardiso library.
//
// FEBio is NOT linked against libpardiso. Instead, the library is located and
// opened with dlopen the first time it is needed. If it is missing, or Panua's
// license check fails, the solver reports itself unavailable and NumCore falls
// back to another default solver (Accelerate on macOS). This lets one binary
// be distributed to users with or without a Panua license.
//
// Search order:
//   1. $PARDISO_LIB_PATH (full path to the library, or the directory holding it)
//   2. ~/Library/Application Support/FEBio/lib   (macOS)
//      ~/.febio/lib                              (Linux)
//   3. ~/panua-pardiso-*/lib  (Panua's default install location; newest first)
//   4. PDL_DEFAULT_PATH  (the PDL_LIB path FEBio was configured with, if any)
//   5. /usr/local/lib, /opt/homebrew/lib (macOS)
//   6. the dynamic loader's default search
#include <dlfcn.h>
#include <glob.h>
#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <mutex>

typedef void (*pardisoinit_fn)(void*, int*, int*, int*, double*, int*);
typedef void (*pardiso_fn)(void*, int*, int*, int*, int*, int*,
                           double*, int*, int*, int*, int*, int*,
                           int*, double*, double*, int*, double*);

static pardisoinit_fn pardisoinit = nullptr;
static pardiso_fn     pardiso     = nullptr;
static std::string    s_pdl_path;     // library that was loaded
static std::string    s_pdl_error;    // reason the library is unavailable
static bool           s_pdl_ok = false;

#ifdef __APPLE__
static const char* PDL_LIBNAME = "libpardiso.dylib";
#else
static const char* PDL_LIBNAME = "libpardiso.so";
#endif

static void load_pardiso_library()
{
	std::vector<std::string> candidates;

	const char* env = getenv("PARDISO_LIB_PATH");
	if (env && *env)
	{
		std::string p(env);
		// accept either a file or a directory
		auto ends_with = [&](const char* sfx) {
			size_t n = strlen(sfx);
			return (p.size() >= n) && (p.compare(p.size() - n, n, sfx) == 0);
		};
		if (!ends_with(".dylib") && !ends_with(".so"))
			p += std::string("/") + PDL_LIBNAME;
		candidates.push_back(p);
	}

	const char* home = getenv("HOME");
	if (home)
	{
#ifdef __APPLE__
		candidates.push_back(std::string(home) + "/Library/Application Support/FEBio/lib/" + PDL_LIBNAME);
#else
		candidates.push_back(std::string(home) + "/.febio/lib/" + PDL_LIBNAME);
#endif

		// Panua's default install location: ~/panua-pardiso-<yyyymmdd>-<platform>/lib.
		// The date sorts lexicographically, so try the newest release first.
		std::string pattern = std::string(home) + "/panua-pardiso-*/lib/" + PDL_LIBNAME;
		glob_t g;
		if (glob(pattern.c_str(), 0, nullptr, &g) == 0)
		{
			std::vector<std::string> found(g.gl_pathv, g.gl_pathv + g.gl_pathc);
			std::sort(found.rbegin(), found.rend());
			candidates.insert(candidates.end(), found.begin(), found.end());
		}
		globfree(&g);
	}

#ifdef PDL_DEFAULT_PATH
	// the PDL_LIB path FEBio was configured with (only useful on the build machine)
	candidates.push_back(PDL_DEFAULT_PATH);
#endif
	candidates.push_back(std::string("/usr/local/lib/") + PDL_LIBNAME);
#ifdef __APPLE__
	candidates.push_back(std::string("/opt/homebrew/lib/") + PDL_LIBNAME);
#endif
	candidates.push_back(PDL_LIBNAME);

	void* h = nullptr;
	std::string lastErr;
	for (const std::string& c : candidates)
	{
		h = dlopen(c.c_str(), RTLD_NOW | RTLD_LOCAL);
		if (h) { s_pdl_path = c; break; }
		const char* e = dlerror();
		// remember the first error from a file that exists but failed to load
		// (e.g. code-signature or quarantine rejection, missing dependency)
		if (e && lastErr.empty() && strstr(e, "no such file") == nullptr && strstr(e, "not found") == nullptr)
			lastErr = e;
	}

	if (h == nullptr)
	{
		s_pdl_error = std::string(PDL_LIBNAME) + " not found";
		if (!lastErr.empty()) s_pdl_error += " (" + lastErr + ")";
		return;
	}

	pardisoinit = (pardisoinit_fn)dlsym(h, "pardisoinit");
	pardiso     = (pardiso_fn)    dlsym(h, "pardiso");
	if ((pardisoinit == nullptr) || (pardiso == nullptr))
	{
		s_pdl_error = s_pdl_path + " does not export pardisoinit/pardiso";
		pardisoinit = nullptr; pardiso = nullptr;
		dlclose(h);
		return;
	}

	// Suppress Panua's license banner, then run the license check.
	setenv("PARDISOLICMESSAGE", "1", 0);
	void* pt[64] = { nullptr };
	int iparm[64] = { 0 };
	double dparm[64] = { 0.0 };
	int mtype = 11, solver = 0, error = 0;
	pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);
	switch (error)
	{
	case 0: s_pdl_ok = true; break;
	case -10: s_pdl_error = "no Panua Pardiso license file found"; break;
	case -11: s_pdl_error = "Panua Pardiso license has expired"; break;
	case -12: s_pdl_error = "Panua Pardiso license: wrong username or hostname"; break;
	default : s_pdl_error = "pardisoinit failed with error " + std::to_string(error); break;
	}
}

static bool pardiso_available()
{
	static std::once_flag once;
	std::call_once(once, load_pardiso_library);
	return s_pdl_ok;
}

#else
/* PARDISO prototype. */
extern "C" void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
extern "C" void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                             double *, int    *,    int *, int *,   int *, int *,
                             int *, double *, double *, int *, double *);
static bool pardiso_available() { return true; }
#endif

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::IsAvailable()
{
	return pardiso_available();
}

std::string PardisoProjectSolver::UnavailableReason()
{
#ifdef PARDISODL_RUNTIME
	pardiso_available();
	return s_pdl_error;
#else
	return std::string();
#endif
}

//-----------------------------------------------------------------------------
// print pardiso error message
void print_err_pdl(int nerror)
{
    switch (-nerror)
    {
        case 1: fprintf(stderr, "Inconsistent input\n"); break;
        case 2: fprintf(stderr, "Not enough memory\n"); break;
        case 3: fprintf(stderr, "Reordering problem\n"); break;
        case 4: fprintf(stderr, "Zero pivot, numerical fact. or iterative refinement problem\n"); break;
        case 5: fprintf(stderr, "Unclassified (internal) error\n"); break;
        case 6: fprintf(stderr, "Preordering failed\n"); break;
        case 7: fprintf(stderr, "Diagonal matrix problem\n"); break;
        case 8: fprintf(stderr, "32-bit integer overflow problem\n"); break;
        default:
            fprintf(stderr, " Unknown\n");
    }
}

//////////////////////////////////////////////////////////////
// PardisoProjectSolver
//////////////////////////////////////////////////////////////

BEGIN_FECORE_CLASS(PardisoProjectSolver, LinearSolver)
ADD_PARAMETER(m_print_cn, "print_condition_number");
ADD_PARAMETER(m_iparm3  , "precondition");
ADD_PARAMETER(m_pivot_perturbation, "pivot_perturbation");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
PardisoProjectSolver::PardisoProjectSolver(FEModel* fem) : LinearSolver(fem), m_pA(0)
{
    m_print_cn = false;
    m_mtype = -2;
    m_iparm3 = false;
    m_pivot_perturbation = 0;  // try 14 when problem does not converge
    m_isFactored = false;
    
    /* If both PARDISO AND PARDISODL are defined, print a warning */
    // Linking both MKL and Panua at build time makes their pardiso/pardisoinit
    // symbols clash. Run-time loading resolves Panua's symbols from its own
    // library handle (dlsym), so both can coexist.
#if defined(PARDISO) && !defined(PARDISODL_RUNTIME)
    fprintf(stderr, "WARNING: The pardiso-project version of the Pardiso solver is being used\n\n");
    exit(1);
#endif
}

//-----------------------------------------------------------------------------
PardisoProjectSolver::~PardisoProjectSolver()
{
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::PrintConditionNumber(bool b)
{
    m_print_cn = b;
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::UseIterativeFactorization(bool b)
{
    m_iparm3 = b;
}

//-----------------------------------------------------------------------------
SparseMatrix* PardisoProjectSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    // allocate the correct matrix format depending on matrix symmetry type
    switch (ntype)
    {
        case REAL_SYMMETRIC     : m_mtype = -2; m_pA = new CompactSymmMatrix(1); break;
        case REAL_UNSYMMETRIC   : m_mtype = 11; m_pA = new CRSSparseMatrix(1); break;
        case REAL_SYMM_STRUCTURE: m_mtype =  1; m_pA = new CRSSparseMatrix(1); break;
        default:
            assert(false);
            m_pA = nullptr;
    }
    
    return m_pA;
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (m_pA && m_isFactored) Destroy();
    m_pA = dynamic_cast<CompactMatrix*>(pA);
    m_mtype = -2;
    if (dynamic_cast<CRSSparseMatrix*>(pA)) m_mtype = 11;
    return (m_pA != nullptr);
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::PreProcess()
{
    m_iparm[0] = 0; /* Use default values for parameters */
//    for (int i=0; i<64; ++i) m_pt[i] = nullptr;
    
    //fprintf(stderr, "In PreProcess\n");
    assert(m_isFactored == false);

    if (!IsAvailable())
    {
        feLogError("pardiso-project solver is unavailable: %s", UnavailableReason().c_str());
        return false;
    }

    int error = 0;
    int solver = 0; /* use sparse direct solver */
    pardisoinit(m_pt, &m_mtype, &solver, m_iparm, m_dparm, &error);
    if (error != 0)
    {
        feLogError("pardisoinit failed with error %d", error);
        return false;
    }

    /* Pivot perturbation threshold 10^-m_pivot_perturbation (iparm(10)).
     The Panua defaults (12 for nonsymmetric, 8 for symmetric) can perturb many
     pivots of badly scaled matrices (e.g. fluid-FSI), degrading the solution.
     A value of 0 keeps the Panua default. */
    if (m_pivot_perturbation > 0) m_iparm[9] = m_pivot_perturbation;
    
    m_n = m_pA->Rows();
    m_nnz = m_pA->NonZeroes();
    m_nrhs = 1;
    
    /* Numbers of processors, value of OMP_NUM_THREADS (if available) or
     total number of processors on machine */
    const auto processor_count = std::thread::hardware_concurrency();
    int num_procs = processor_count;
    char* var = getenv("OMP_NUM_THREADS");
    if (var != NULL)
        sscanf( var, "%d", &num_procs );
    feLog("Number of processors: %d\n",num_procs);
    if (num_procs == 0) {
        feLogError("Specify number of processors in environmental variable OMP_NUM_THREADS!");
        return false;
    }
    m_iparm[2]  = num_procs;

    m_maxfct = 1;    /* Maximum number of numerical factorizations */
    m_mnum = 1;    /* Which factorization to use */
    
    m_msglvl = 0;    /* 0 Suppress printing, 1 Print statistical information */
    
    return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::Factor()
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    
    // ------------------------------------------------------------------------------
    // Reordering and Symbolic Factorization.  This step also allocates all memory
    // that is necessary for the factorization.
    // ------------------------------------------------------------------------------
    
    int phase = 11;
    
    int error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during symbolic factorization: ");
        print_err_pdl(error);
        exit(2);
    }
    
    // ------------------------------------------------------------------------------
    // This step does the factorization
    // ------------------------------------------------------------------------------
    
    phase = 22;
    
    m_iparm[3] = (m_iparm3 ? 61 : 0);
    error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during factorization: ");
        print_err_pdl(error);
        return false;
    }
    
    // calculate and print the condition number
    if (m_print_cn)
    {
        double c = condition_number();
        feLog("\tcondition number (est.) ................... : %lg\n\n", c);
    }
    
    m_isFactored = true;
    
    return true;
}

//-----------------------------------------------------------------------------
bool PardisoProjectSolver::BackSolve(double* x, double* b)
{
    // make sure we have work to do
    if (m_pA->Rows() == 0) return true;
    
    int phase = 33;
    
    m_iparm[7] = 1;    /* Maximum number of iterative refinement steps */
    
    int error = 0;
    pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, m_pA->Values(), m_pA->Pointers(), m_pA->Indices(),
            NULL, &m_nrhs, m_iparm, &m_msglvl, b, x, &error, m_dparm);
    
    if (error)
    {
        fprintf(stderr, "\nERROR during solution: ");
        print_err_pdl(error);
        exit(3);
    }
    
    // update stats
    UpdateStats(1);
    
    return true;
}

//-----------------------------------------------------------------------------
// This algorithm (naively) estimates the condition number. It is based on the observation that
// for a linear system of equations A.x = b, the following holds
// || A^-1 || >= ||x||.||b||
// Thus the condition number can be estimated by
// c = ||A||.||A^-1|| >= ||A|| . ||x|| / ||b||
// This algorithm tries for some random b vectors with norm ||b||=1 to maxize the ||x||.
// The returned value will be an underestimate of the condition number
double PardisoProjectSolver::condition_number()
{
    // This assumes that the factorization is already done!
    int N = m_pA->Rows();
    
    // get the norm of the matrix
    double normA = m_pA->infNorm();
    
    // estimate the norm of the inverse of A
    double normAi = 0.0;
    
    // choose max iterations
    int iters = (N < 50 ? N : 50);
    
    vector<double> b(N, 0), x(N, 0);
    for (int i = 0; i < iters; ++i)
    {
        // create a random vector
        NumCore::randomVector(b, -1.0, 1.0);
        for (int j = 0; j < N; ++j) b[j] = (b[j] >= 0.0 ? 1.0 : -1.0);
        
        // calculate solution
        BackSolve(&x[0], &b[0]);
        
        double normb = NumCore::infNorm(b);
        double normx = NumCore::infNorm(x);
        if (normx > normAi) normAi = normx;
        
        int pct = (100 * i) / (iters - 1);
        fprintf(stderr, "calculating condition number: %d%%\r", pct);
    }
    
    double c = normA*normAi;
    return c;
}

//-----------------------------------------------------------------------------
void PardisoProjectSolver::Destroy()
{
    int phase = -1;
    
    int error = 0;
    
    if (m_pA && m_pA->Pointers() && m_isFactored)
    {
        pardiso(m_pt, &m_maxfct, &m_mnum, &m_mtype, &phase, &m_n, NULL, m_pA->Pointers(), m_pA->Indices(),
                NULL, &m_nrhs, m_iparm, &m_msglvl, NULL, NULL, &error, m_dparm);
    }
    m_isFactored = false;
}

#else
bool PardisoProjectSolver::IsAvailable() { return false; }
std::string PardisoProjectSolver::UnavailableReason() { return "FEBio was built without pardiso-project support"; }

BEGIN_FECORE_CLASS(PardisoProjectSolver, LinearSolver)
	ADD_PARAMETER(m_print_cn, "print_condition_number");
	ADD_PARAMETER(m_iparm3, "precondition");
	ADD_PARAMETER(m_pivot_perturbation, "pivot_perturbation");
END_FECORE_CLASS();

PardisoProjectSolver::PardisoProjectSolver(FEModel* fem) : LinearSolver(fem) { m_pivot_perturbation = 14; }
PardisoProjectSolver::~PardisoProjectSolver() {}
bool PardisoProjectSolver::PreProcess() { return false; }
bool PardisoProjectSolver::Factor() { return false; }
bool PardisoProjectSolver::BackSolve(double* x, double* y) { return false; }
void PardisoProjectSolver::Destroy() {}
SparseMatrix* PardisoProjectSolver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool PardisoProjectSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
void PardisoProjectSolver::PrintConditionNumber(bool b) {}
double PardisoProjectSolver::condition_number() { return 0; }
void PardisoProjectSolver::UseIterativeFactorization(bool b) {}
#endif
