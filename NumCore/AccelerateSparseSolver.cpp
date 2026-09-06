/*This file is part of the FEBio source code and is licensed under the MIT license
 listed below.

 See Copyright-FEBio.txt for details.

 Copyright (c) 2022 University of Utah, The Trustees of Columbia University in
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
#include "AccelerateSparseSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>
#include <cmath>
#include <vector>
#include <string>
#include <string.h>

#ifdef WIN32
#define fecmpi _strnicmp
#else
#define fecmpi strncasecmp
#endif

// Sparse LU was added to Accelerate in macOS 15.5 / iOS 18.5. The
// SparseFactorizationLU* constants do not exist in older SDKs, so guard their
// use at compile time as well as at run time.
#if (defined(__MAC_OS_X_VERSION_MAX_ALLOWED) && (__MAC_OS_X_VERSION_MAX_ALLOWED >= 150500)) \
 || (defined(__IPHONE_OS_VERSION_MAX_ALLOWED) && (__IPHONE_OS_VERSION_MAX_ALLOWED >= 180500))
#define FEBIO_ACCELERATE_HAS_LU 1
#else
#define FEBIO_ACCELERATE_HAS_LU 0
#endif

#ifdef HAS_ACCEL

class AccelerateSparseSolver::Implementation
{
public:
    FEModel*            m_fem;

    CompactMatrix*      m_pA;
    int                 m_mtype; // matrix type: -2 = symmetric, 11 = unsymmetric

#ifdef __APPLE__
    SparseMatrixStructure SMS;
    SparseMatrix_Double A;
    SparseOpaqueSymbolicFactorization ASS;
    SparseSymbolicFactorOptions SSFO;
    SparseNumericFactorOptions SNFO;
    SparseOpaqueFactorization_Double ASF;
    long* colS;
    std::vector<long> pointers;
    SparseIterativeStatus_t SStatus;
#endif

    // partition-aware diagonal scaling
    std::vector<double> m_scale;         // D, one entry per equation
    std::vector<double> m_scaledValues;  // values of A' = D A D
    std::vector<double> m_bs;            // scratch for the scaled right-hand side

    // Matrix data
    int m_n, m_nnz, m_nrhs;

    bool    m_print_cn;     // estimate and print the condition number
    bool    m_iparm3;       // use iterative method
    bool    m_isSymbolic;   // symbolic factorization exists (must be cleaned up)
    bool    m_isFactored;   // numeric factorization exists (must be cleaned up)
    int     m_ftype;        // factorization type
    int     m_ordmthd;      // order method
    int     m_itrmthd;      // iterative method
    int     m_precond;      // preconditioner for the iterative methods
    int     m_maxiter;      // max number of iterations
    int     m_nvec;         // nr. of orthogonal vectors kept by GMRES (0 = Accelerate default of 16)
    double  m_rtol;         // relative tolerance
    double  m_atol;         // absolute tolerance
    int     m_print_level;  // print level for status updates
    int     m_scaling;      // partition scaling (PartitionScalingOff / PartitionScalingOn)
    int     m_strategy;     // SSDirect / SSIterative / SSFactorPreconditioned
    int     m_pcMaxAge;     // refactor after this many preconditioned solves (0 = never on age alone)
    int     m_pcAge;        // solves performed with the current factorization
    bool    m_pcRefactor;   // force a refactorization on the next Factor()

    // String forms of the option parameters, as they appear in the input file.
    //
    // These are read as plain strings rather than as FEBio enum-ints because
    // the enum machinery is not reachable on every parameter-reading path (the
    // febio.xml config reader in particular), which left symbolic values such
    // as "qr" silently unusable. Reading a string and resolving it ourselves in
    // ResolveParameters() works identically no matter which reader runs, and
    // still accepts the plain integers that older input files use.
    //
    // An empty string means "not specified" - the corresponding integer member
    // keeps whatever the constructor or a programmatic setter put there.
    std::string m_sftype;
    std::string m_sordmthd;
    std::string m_sitrmthd;
    std::string m_sprecond;
    std::string m_sstrategy;
    std::string m_sscaling;
    std::string m_sfscaling;

    int     m_fscaling;      // Accelerate's own scaling during numeric factorization
    double  m_pivotTol;      // pivot tolerance for threshold partial pivoting, [0, 0.5]

public:
    //! Resolved solve strategy. The legacy boolean "iterative" parameter is
    //! folded in here so that old input files keep working.
    int strategy() const
    {
        if (m_iparm3 && (m_strategy == SSDirect)) return SSIterative;
        return m_strategy;
    }

public:
    Implementation()
    {
        m_print_cn    = false;
        m_mtype       = -2;
        m_iparm3      = false;
        m_isSymbolic  = false;
        m_isFactored  = false;
        m_ftype       = FTSparseFactorizationDefault;   // == 7, a valid enum index
        m_ordmthd     = OMSparseOrderDefault;           // == 3, a valid enum index
        m_itrmthd     = ITSparseLSMR;
        m_precond     = PCDiagonal;
        m_pA          = nullptr;
        m_maxiter     = 1000;
        m_nvec        = 0;
        m_rtol        = 1e-8;   // NOTE: must be double. Was int, which truncated to 0.
        m_atol        = 1e-12;  // NOTE: must be double. Was int, which truncated to 0.
        m_print_level = 0;
        m_scaling     = PartitionScalingOff;
        m_fscaling    = FSDefault;
        m_pivotTol    = -1.0;          // < 0 means "leave Accelerate's default alone"
        m_strategy    = SSDirect;
        m_pcMaxAge    = 0;
        m_pcAge       = 0;
        m_pcRefactor  = false;
        m_n = m_nnz = m_nrhs = 0;
        colS = nullptr;
    }
};

//////////////////////////////////////////////////////////////
// AccelerateSparseSolver
//////////////////////////////////////////////////////////////

// The option parameters are registered as STRINGS, not as FEBio enum-ints.
//
// FEBio does have enum support (FEParam::setEnums plus fexml::enumValue), and
// it works on the model-file path. But it is not reachable on every path that
// reads a parameter list -- notably the febio.xml config reader -- so symbolic
// values such as "qr" were being ignored there. FE_PARAM_STD_STRING is handled
// identically by every reader, so resolving the names ourselves in
// ResolveParameters() makes them work everywhere.
//
// ResolveParameters() also accepts the bare integers that older input files
// use, so nothing that worked before stops working.
BEGIN_FECORE_CLASS(AccelerateSparseSolver, LinearSolver)
    ADD_PARAMETER(imp->m_print_cn   , "print_condition_number");
    ADD_PARAMETER(imp->m_iparm3     , "iterative");   // legacy; folded into "strategy"
    ADD_PARAMETER(imp->m_sstrategy  , "strategy");
    ADD_PARAMETER(imp->m_pcMaxAge   , "preconditioner_max_age");
    ADD_PARAMETER(imp->m_sftype     , "factorization");
    ADD_PARAMETER(imp->m_sordmthd   , "order_method");
    ADD_PARAMETER(imp->m_sitrmthd   , "iterative_method");
    ADD_PARAMETER(imp->m_sprecond   , "preconditioner");
    ADD_PARAMETER(imp->m_maxiter    , "max_iter");
    ADD_PARAMETER(imp->m_nvec       , "nvec");
    ADD_PARAMETER(imp->m_print_level, "print_level");
    ADD_PARAMETER(imp->m_rtol       , "rtol");
    ADD_PARAMETER(imp->m_atol       , "atol");
    ADD_PARAMETER(imp->m_sscaling   , "partition_scaling");
    ADD_PARAMETER(imp->m_sfscaling  , "factor_scaling");
    ADD_PARAMETER(imp->m_pivotTol   , "pivot_tolerance");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
// Name -> value tables for the string-valued options.
namespace {

struct EnumEntry { const char* name; int value; };

const EnumEntry g_ftype[] = {
    { "cholesky",       AccelerateSparseSolver::FTSparseFactorizationCholesky      },
    { "ldlt",           AccelerateSparseSolver::FTSparseFactorizationLDLT          },
    { "ldlt_unpivoted", AccelerateSparseSolver::FTSparseFactorizationLDLTUnpivoted },
    { "ldlt_sbk",       AccelerateSparseSolver::FTSparseFactorizationLDLTSBK       },
    { "ldlt_tpp",       AccelerateSparseSolver::FTSparseFactorizationLDLTTPP       },
    { "qr",             AccelerateSparseSolver::FTSparseFactorizationQR            },
    { "cholesky_ata",   AccelerateSparseSolver::FTSparseFactorizationCholeskyAtA   },
    { "lu",             AccelerateSparseSolver::FTSparseFactorizationLU            },
    { "lu_unpivoted",   AccelerateSparseSolver::FTSparseFactorizationLUUnpivoted   },
    { "lu_spp",         AccelerateSparseSolver::FTSparseFactorizationLUSPP         },
    { "lu_tpp",         AccelerateSparseSolver::FTSparseFactorizationLUTPP         },
    { "default",        AccelerateSparseSolver::FTSparseFactorizationDefault       },
    { nullptr, 0 }
};

const EnumEntry g_fscaling[] = {
    { "default",            AccelerateSparseSolver::FSDefault           },
    { "none",               AccelerateSparseSolver::FSNone              },
    { "equilibriation_inf", AccelerateSparseSolver::FSEquilibriationInf },
    { nullptr, 0 }
};

const EnumEntry g_order[] = {
    { "amd",     AccelerateSparseSolver::OMSparseOrderAMD     },
    { "metis",   AccelerateSparseSolver::OMSparseOrderMetis   },
    { "colamd",  AccelerateSparseSolver::OMSparseOrderCOLAMD  },
    { "default", AccelerateSparseSolver::OMSparseOrderDefault },
    { nullptr, 0 }
};

const EnumEntry g_itmethod[] = {
    { "cg",      AccelerateSparseSolver::ITSparseConjugateGradient },
    { "gmres",   AccelerateSparseSolver::ITSparseGMRES             },
    { "dqgmres", AccelerateSparseSolver::ITSparseDQGMRES           },
    { "fgmres",  AccelerateSparseSolver::ITSparseFGMRES            },
    { "lsmr",    AccelerateSparseSolver::ITSparseLSMR              },
    { nullptr, 0 }
};

const EnumEntry g_precond[] = {
    { "none",         AccelerateSparseSolver::PCNone        },
    { "diagonal",     AccelerateSparseSolver::PCDiagonal    },
    { "diag_scaling", AccelerateSparseSolver::PCDiagScaling },
    { nullptr, 0 }
};

const EnumEntry g_strategy[] = {
    { "direct",                AccelerateSparseSolver::SSDirect               },
    { "iterative",             AccelerateSparseSolver::SSIterative            },
    { "factor_preconditioned", AccelerateSparseSolver::SSFactorPreconditioned },
    { nullptr, 0 }
};

const EnumEntry g_scaling[] = {
    { "off", AccelerateSparseSolver::PartitionScalingOff },
    { "on",  AccelerateSparseSolver::PartitionScalingOn  },
    { "0",   AccelerateSparseSolver::PartitionScalingOff },
    { "1",   AccelerateSparseSolver::PartitionScalingOn  },
    { nullptr, 0 }
};

std::string enumList(const EnumEntry* tbl)
{
    std::string s;
    for (const EnumEntry* e = tbl; e->name; ++e) {
        if (!s.empty()) s += ", ";
        s += e->name;
    }
    return s;
}

// Resolve one option. Returns false only when the string is non-empty and
// matches neither a name nor an integer. A blank string leaves 'out' alone.
bool resolveOne(const std::string& sin, const EnumEntry* tbl, const char* paramName, int& out)
{
    const char* ws = " \t\r\n";
    size_t b = sin.find_first_not_of(ws);
    if (b == std::string::npos) return true;            // not specified
    size_t e = sin.find_last_not_of(ws);
    std::string s = sin.substr(b, e - b + 1);

    for (const EnumEntry* p = tbl; p->name; ++p) {
        if ((s.size() == strlen(p->name)) && (fecmpi(s.c_str(), p->name, s.size()) == 0)) {
            out = p->value;
            return true;
        }
    }

    // fall back to a plain integer, so old numeric input files still work
    char* endp = nullptr;
    long v = strtol(s.c_str(), &endp, 10);
    if (endp && (*endp == 0)) { out = (int) v; return true; }

    fprintf(stderr, "\nERROR: unrecognized value \"%s\" for %s.\n", s.c_str(), paramName);
    fprintf(stderr, "Valid values are: %s (or the equivalent integer).\n", enumList(tbl).c_str());
    return false;
}

} // namespace

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::ResolveParameters()
{
    if (!resolveOne(imp->m_sstrategy, g_strategy,  "strategy",          imp->m_strategy)) return false;
    if (!resolveOne(imp->m_sftype,    g_ftype,     "factorization",     imp->m_ftype   )) return false;
    if (!resolveOne(imp->m_sordmthd,  g_order,     "order_method",      imp->m_ordmthd )) return false;
    if (!resolveOne(imp->m_sitrmthd,  g_itmethod,  "iterative_method",  imp->m_itrmthd )) return false;
    if (!resolveOne(imp->m_sprecond,  g_precond,   "preconditioner",    imp->m_precond )) return false;
    if (!resolveOne(imp->m_sscaling,  g_scaling,   "partition_scaling", imp->m_scaling )) return false;
    if (!resolveOne(imp->m_sfscaling, g_fscaling,  "factor_scaling",    imp->m_fscaling)) return false;
    return true;
}

//-----------------------------------------------------------------------------
// Translate the FEBio-facing factorization id into Accelerate's constant.
// Keeping this in one place is what lets the two numbering schemes differ
// safely (FEBio needs 0..n-1 for named enums; Accelerate uses QR = 40,
// CholeskyAtA = 41).
static bool luAvailable()
{
#if FEBIO_ACCELERATE_HAS_LU
    if (__builtin_available(macOS 15.5, iOS 18.5, *)) return true;
#endif
    return false;
}

static bool toSparseFactorization(int ftype, bool bsymm, SparseFactorization_t& out)
{
    if (ftype == AccelerateSparseSolver::FTSparseFactorizationDefault) {
        if (bsymm) ftype = AccelerateSparseSolver::FTSparseFactorizationLDLTTPP;
        else {
            // LU is the right default for a general matrix: for a square system
            // it costs roughly half the flops of QR and produces far less fill,
            // since QR's R factor carries the sparsity of chol(A^T A). Fall back
            // to QR on systems older than macOS 15.5.
            ftype = luAvailable() ? AccelerateSparseSolver::FTSparseFactorizationLU
                                  : AccelerateSparseSolver::FTSparseFactorizationQR;
        }
    }

#if FEBIO_ACCELERATE_HAS_LU
    if (__builtin_available(macOS 15.5, iOS 18.5, *)) {
        switch (ftype) {
            case AccelerateSparseSolver::FTSparseFactorizationLU:
                out = SparseFactorizationLU;           return !bsymm;
            case AccelerateSparseSolver::FTSparseFactorizationLUUnpivoted:
                out = SparseFactorizationLUUnpivoted;  return !bsymm;
            case AccelerateSparseSolver::FTSparseFactorizationLUSPP:
                out = SparseFactorizationLUSPP;        return !bsymm;
            case AccelerateSparseSolver::FTSparseFactorizationLUTPP:
                out = SparseFactorizationLUTPP;        return !bsymm;
            default: break;
        }
    }
    else
#endif
    {
        switch (ftype) {
            case AccelerateSparseSolver::FTSparseFactorizationLU:
            case AccelerateSparseSolver::FTSparseFactorizationLUUnpivoted:
            case AccelerateSparseSolver::FTSparseFactorizationLUSPP:
            case AccelerateSparseSolver::FTSparseFactorizationLUTPP:
                fprintf(stderr, "\nSparse LU requires macOS 15.5 / iOS 18.5 or later"
                                " (both to build and to run). Use qr instead.\n");
                return false;
            default: break;
        }
    }

    switch (ftype) {
        case AccelerateSparseSolver::FTSparseFactorizationCholesky:
            out = SparseFactorizationCholesky;        return bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationLDLT:
            out = SparseFactorizationLDLT;            return bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationLDLTUnpivoted:
            out = SparseFactorizationLDLTUnpivoted;   return bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationLDLTSBK:
            out = SparseFactorizationLDLTSBK;         return bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationLDLTTPP:
            out = SparseFactorizationLDLTTPP;         return bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationQR:
            out = SparseFactorizationQR;              return !bsymm;
        case AccelerateSparseSolver::FTSparseFactorizationCholeskyAtA:
            out = SparseFactorizationCholeskyAtA;     return !bsymm;
        default:
            return false;
    }
}

//-----------------------------------------------------------------------------
AccelerateSparseSolver::AccelerateSparseSolver(FEModel* fem) : LinearSolver(fem), imp(new AccelerateSparseSolver::Implementation)
{
    imp->m_fem = fem;
}

//-----------------------------------------------------------------------------
AccelerateSparseSolver::~AccelerateSparseSolver()
{
    Destroy();
    delete imp;         // was leaked
    imp = nullptr;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::PrintConditionNumber(bool b)
{
    imp->m_print_cn = b;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::UseIterativeFactorization(bool b)
{
    imp->m_iparm3 = b;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetFactorizationType(int ftype)
{
    imp->m_ftype = ftype;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetOrderMethod(int order)
{
    imp->m_ordmthd = order;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetIterativeMethod(int method)
{
    imp->m_itrmthd = method;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetPartitionScaling(PartitionScaling s)
{
    imp->m_scaling = s;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetSolveStrategy(int strategy)
{
    imp->m_strategy = strategy;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::SetPrintLevel(int printlevel)
{
    imp->m_print_level = printlevel;
}

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::IsIterative() const
{
    // factor_preconditioned still solves iteratively, it just happens to carry
    // a direct factorization around as the preconditioner.
    return (imp->strategy() != SSDirect);
}

//-----------------------------------------------------------------------------
// Preconditioner callback: Y = P X with P the stored factorization of A.
//
// Accelerate documents this as computing Y = P X where P approximates A^-1, so
// "applying" the preconditioner means doing a triangular solve with the stored
// factors. The factorization is deliberately allowed to be out of date.
void AccelerateSparseSolver::ApplyStoredFactor(void* mem, enum CBLAS_TRANSPOSE trans,
                                               DenseMatrix_Double X, DenseMatrix_Double Y)
{
    Implementation* im = (Implementation*) mem;

    // Fall back to the identity preconditioner when we cannot apply the stored
    // factors. GMRES stays mathematically correct with P = I; it just converges
    // slowly, which is far better than returning garbage.
    //
    // The transpose case is one of these: SparseSolve() on a stored
    // factorization solves A x = b only. A transpose solve would have to be
    // assembled by hand from the subfactors (Q/R/P for QR, P/Q/S_r/S_c for LU).
    // GMRES never asks for it (LSMR does), so this branch should not fire for
    // the intended configuration.
    const bool bok = im->m_isFactored && (trans == CblasNoTrans);

    if (!bok) {
        for (int c = 0; c < X.columnCount; ++c)
            for (int r = 0; r < X.rowCount; ++r)
                Y.data[(size_t)c*Y.columnStride + r] = X.data[(size_t)c*X.columnStride + r];
        return;
    }

    // SparseSolve(F, B, X) solves A X = B, so the RHS is our input X and the
    // result lands in our output Y.
    SparseSolve(im->ASF, X, Y);
}

//-----------------------------------------------------------------------------
// NOTE: this callback is installed both on SSFO.reportError (symbolic
// factorization) and on the iterative option structs, so the message must not
// claim to come from either one. It previously said "ERROR during iterative
// solve", which made symbolic-factorization failures look like iterative-solver
// problems.
void AccelerateSparseSolver::MyReportError(const char* message) {
    fprintf(stderr, "\n\tAccelerate sparse solver error: %s\n", message);
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::MyReportStatus(const char* message) {
    fprintf(stdout, "iterative solver status: %s",message);
    return;
}

//-----------------------------------------------------------------------------
SparseMatrix* AccelerateSparseSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    // Allocate the correct matrix format depending on matrix symmetry type.
    //
    // NOTE: REAL_SYMMETRIC must use CompactSymmMatrix (column-major, LOWER
    // triangle only). The previous code allocated a CCSSparseMatrix here, which
    // stores BOTH triangles, and then told Accelerate the matrix was
    // SparseSymmetric. Accelerate transposes and SUMS any entry it finds in the
    // wrong triangle, so every off-diagonal was silently doubled. The
    // factorization succeeded and returned a plausible but wrong answer.
    //
    // REAL_SYMM_STRUCTURE is only structurally symmetric; numerically it is a
    // general matrix and must be stored and factored as such.
    switch (ntype)
    {
        case REAL_SYMMETRIC     : imp->m_mtype = -2; imp->m_pA = new CompactSymmMatrix(0); break;
        case REAL_UNSYMMETRIC   : imp->m_mtype = 11; imp->m_pA = new CCSSparseMatrix(0);   break;
        case REAL_SYMM_STRUCTURE: imp->m_mtype = 11; imp->m_pA = new CCSSparseMatrix(0);   break;
        default:
            assert(false);
            imp->m_pA = nullptr;
    }

    return imp->m_pA;
}

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (imp->m_pA && (imp->m_isFactored || imp->m_isSymbolic)) Destroy();

    imp->m_pA = dynamic_cast<CompactMatrix*>(pA);
    if (imp->m_pA == nullptr) return false;

    // Derive the matrix type from the storage itself.
    //
    // The previous version did:
    //     m_mtype = -2;
    //     if (dynamic_cast<CCSSparseMatrix*>(pA)) m_mtype = 11;
    // but CreateSparseMatrix() handed out a CCSSparseMatrix for *every* matrix
    // type, so this branch always fired and every matrix was classified as
    // unsymmetric -- discarding the symmetry information that
    // CreateSparseMatrix() had correctly established. Asking for a Cholesky or
    // LDLT factorization then produced
    //     "Incorrect factorization type for non-symmetric matrix!"
    // on a matrix that really was symmetric.
    imp->m_mtype = imp->m_pA->isSymmetric() ? -2 : 11;

    return true;
}

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::PreProcess()
{
    assert(imp->m_isFactored == false);

    if (imp->m_pA == nullptr) return false;

    // Turn the string options ("qr", "gmres", "on", ...) into integers before
    // anything reads them.
    if (!ResolveParameters()) return false;

    if (imp->m_print_level > 0) {
        feLog("\tAccelerate solver: strategy = %d, factorization = %d, order = %d,\n"
              "\t                   iterative method = %d, preconditioner = %d, scaling = %d\n",
              imp->strategy(), imp->m_ftype, imp->m_ordmthd,
              imp->m_itrmthd, imp->m_precond, imp->m_scaling);
    }

    imp->m_n    = imp->m_pA->Rows();
    imp->m_nnz  = imp->m_pA->NonZeroes();
    imp->m_nrhs = 1;

    const int ncols = imp->m_pA->Columns();

    // Accelerate expects 0-based indices.
    if (imp->m_pA->Offset() != 0) {
        fprintf(stderr, "\nERROR during preprocessing: ");
        fprintf(stderr, "\nAccelerate solver requires a 0-based matrix (offset = %d).\n", imp->m_pA->Offset());
        return false;
    }

    // Accelerate's SparseMatrixStructure is column-based. Handing it the row
    // pointers of a CRS matrix would silently factor the transpose.
    if (imp->m_pA->isRowBased()) {
        fprintf(stderr, "\nERROR during preprocessing: ");
        fprintf(stderr, "\nAccelerate solver requires a column-based (CCS) matrix.\n");
        return false;
    }

    // columnStarts has columnCount+1 entries -- NOT nnz.
    //
    // The previous code did pointers.resize(m_nnz) and then copied m_nnz values
    // out of Pointers(), which only has ncols+1 elements. For any realistic
    // matrix nnz >> ncols+1, so this read far past the end of the pointer array.
    imp->pointers.resize(ncols + 1);
    for (int i = 0; i <= ncols; ++i) imp->pointers[i] = imp->m_pA->Pointers()[i];
    imp->colS = &imp->pointers[0];

    if (__builtin_available(macOS 10.13, *)) {
        imp->SMS.columnCount  = ncols;
        imp->SMS.rowCount     = imp->m_pA->Rows();
        imp->SMS.rowIndices   = imp->m_pA->Indices();
        imp->SMS.columnStarts = imp->colS;
        imp->SMS.blockSize    = 1;

        // Zero the whole attribute bitfield before touching individual members.
        // SparseAttributes_t contains _reserved (documented as "must be zero")
        // and _allocatedBySparse, which tells SparseCleanup() whether it owns
        // the storage. Leaving those bits as stack garbage risks Accelerate
        // trying to free memory that belongs to FEBio.
        imp->SMS.attributes = SparseAttributes_t{};

        // Both SSDirect and SSFactorPreconditioned need a symbolic
        // factorization; only SSIterative skips it entirely.
        if (imp->strategy() != SSIterative) {
            imp->SSFO = SparseSymbolicFactorOptions{};
            imp->SSFO.control              = SparseDefaultControl;
            imp->SSFO.order                = nullptr;
            imp->SSFO.ignoreRowsAndColumns = nullptr;
            imp->SSFO.malloc               = malloc;
            imp->SSFO.free                 = free;
            imp->SSFO.reportError          = MyReportError;

            imp->SNFO = SparseNumericFactorOptions{};
            imp->SNFO.control = SparseDefaultControl;
            imp->SNFO.scaling = nullptr;

            // Accelerate's own equilibration. Note SparseScalingUser with a
            // null scaling array means "no scaling at all", which is how FSNone
            // is expressed. For LU the framework default is already no scaling,
            // so equilibriation_inf is the interesting setting here.
            switch (imp->m_fscaling) {
                case FSNone:              imp->SNFO.scalingMethod = SparseScalingUser;              break;
                case FSEquilibriationInf: imp->SNFO.scalingMethod = SparseScalingEquilibriationInf; break;
                case FSDefault:
                default:                  imp->SNFO.scalingMethod = SparseScalingDefault;           break;
            }

            // Pivot tolerance for the threshold-partial-pivoting variants
            // (LDLTTPP, LUTPP, and the LU/LDLT defaults which currently are
            // TPP). Accelerate clamps this to [0, 0.5]; a negative value here
            // means the user did not set it, so leave the struct's default.
            if (imp->m_pivotTol >= 0.0) imp->SNFO.pivotTolerance = imp->m_pivotTol;

            const bool bsymm = (imp->m_mtype != 11);

            // Matrix attributes.
            if (bsymm) {
                // CompactSymmMatrix stores the lower triangle in column-major
                // order, so tell Accelerate exactly that.
                imp->SMS.attributes.kind     = SparseSymmetric;
                imp->SMS.attributes.triangle = SparseLowerTriangle;
            }
            else {
                imp->SMS.attributes.kind = SparseOrdinary;   // 'triangle' is ignored
            }

            // Resolve the factorization BEFORE the ordering: which orderings are
            // legal depends on which factorization was chosen.
            SparseFactorization_t ftype;
            if (!toSparseFactorization(imp->m_ftype, bsymm, ftype)) {
                fprintf(stderr, "\nIncorrect factorization type (%d) for %s matrix!",
                        imp->m_ftype, bsymm ? "symmetric" : "non-symmetric");
                if (bsymm)
                    fprintf(stderr, "\nValid choices are: cholesky, ldlt, ldlt_unpivoted, ldlt_sbk, ldlt_tpp, default.\n");
                else {
                    fprintf(stderr, "\nValid choices are: lu, lu_spp, lu_tpp, lu_unpivoted, qr, cholesky_ata, default.");
                    fprintf(stderr, "\n(lu* requires macOS 15.5 / iOS 18.5 or later.)");
                    fprintf(stderr, "\ncholesky and the ldlt* variants require a symmetric matrix.\n");
                }
                return false;
            }

            // COLAMD computes an ordering for A^T A. That is exactly what the
            // QR-based factorizations want, and it is rejected by everything
            // else: by the symmetric factorizations (documented) and by LU,
            // which needs an ordering of A itself (Dulmage-Mendelsohn, with
            // separate row and column permutations). Requesting it for LU gives
            //
            //   The specified options.orderMethod (4) is not supported for
            //   this factorization type (83)
            //
            // where 83 is SparseFactorizationLUTPP, the type that "lu" resolves
            // to since the default LU is currently TPP.
            const bool bQRbased = (ftype == SparseFactorizationQR)
                               || (ftype == SparseFactorizationCholeskyAtA);

            switch (imp->m_ordmthd) {
                case OMSparseOrderDefault:
                    // Let Accelerate choose. It knows what suits the resolved
                    // factorization type, and hard-coding a guess here is what
                    // produced the COLAMD/LU mismatch in the first place. (The
                    // header's claim that SparseOrderDefault means "COLAMD for
                    // unsymmetric factorizations" predates LU.)
                    imp->SSFO.orderMethod = SparseOrderDefault;
                    break;

                case OMSparseOrderAMD:
                    imp->SSFO.orderMethod = SparseOrderAMD;
                    break;

                case OMSparseOrderMetis:
                    imp->SSFO.orderMethod = SparseOrderMetis;
                    break;

                case OMSparseOrderCOLAMD:
                    if (!bQRbased) {
                        fprintf(stderr, "\nOrder method 'colamd' is only valid for the qr and cholesky_ata"
                                        "\nfactorizations (it orders A^T A). Use amd, metis, or default.\n");
                        return false;
                    }
                    imp->SSFO.orderMethod = SparseOrderCOLAMD;
                    break;

                default:
                    fprintf(stderr, "\nIncorrect order method (%d)!", imp->m_ordmthd);
                    fprintf(stderr, "\nValid choices are: amd, metis, colamd, default.\n");
                    return false;
            }

            imp->ASS = SparseFactor(ftype, imp->SMS, imp->SSFO);

            if (imp->ASS.status != SparseStatusOK) {
                fprintf(stderr, "\n\tERROR during symbolic factorization (status = %d)\n", (int)imp->ASS.status);
                return false;
            }
            imp->m_isSymbolic = true;
            imp->m_pcRefactor = true;   // nothing factored yet
            imp->m_pcAge      = 0;
        }
        else {
            // For iterative solving we don't factorize, but the matrix must
            // still be described as a general matrix so that the matrix-vector
            // products use both triangles.
            imp->SMS.attributes.kind = SparseOrdinary;
            if (imp->m_mtype != 11) {
                imp->SMS.attributes.kind     = SparseSymmetric;
                imp->SMS.attributes.triangle = SparseLowerTriangle;
            }
            return true;
        }
    }
    else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during preprocessing: ");
        fprintf(stderr, "\naccelerate solver not available for macOS earlier than 10.13!");
        return false;
    }

    return LinearSolver::PreProcess();
}

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::Factor()
{
    // make sure we have work to do
    if (imp->m_pA->Rows() == 0) return true;

    // ------------------------------------------------------------------------------
    // This step does the factorization
    // ------------------------------------------------------------------------------

    if (__builtin_available(macOS 10.13, *)) {

        // Compute D and the scaled values A' = D A D. When scaling is off this
        // fills m_scale with 1.0 and leaves m_scaledValues empty.
        if (!ComputeAndApplyScaling()) return false;

        // create the sparse matrix
        //
        // Point Accelerate at the SCALED values. The previous version computed
        // the scaling and then handed over m_pA->Values(), so the scaling was
        // discarded here and applied nowhere.
        const bool bscale = (imp->m_scaling == PartitionScalingOn);
        imp->A.data      = bscale ? imp->m_scaledValues.data() : imp->m_pA->Values();
        imp->A.structure = imp->SMS;

        const int strat = imp->strategy();

        if (strat == SSIterative) {
            // no factorization at all
            return true;
        }
        else if (strat == SSFactorPreconditioned) {
            // Refactor only when the stored preconditioner is missing, has been
            // explicitly invalidated (a previous solve failed to converge), or
            // has been used more than preconditioner_max_age times. Otherwise
            // keep the stale factors -- that is the entire point of this mode.
            const bool bneed = (!imp->m_isFactored)
                            || imp->m_pcRefactor
                            || ((imp->m_pcMaxAge > 0) && (imp->m_pcAge >= imp->m_pcMaxAge));

            if (bneed) {
                if (!FactorNow()) return false;
                if (imp->m_print_level > 0)
                    feLog("\trefactored preconditioner (age was %d)\n", imp->m_pcAge);
                imp->m_pcAge     = 0;
                imp->m_pcRefactor = false;
            }
            return true;
        }
        else {
            // SSDirect: factor every time
            if (!FactorNow()) return false;
        }
    }
    else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during factorization: ");
        fprintf(stderr, "\nAccelerate solver not available for macOS earlier than 10.13!");
        return false;
    }

    // calculate and print the condition number
    // (m_isFactored is set inside FactorNow(), before this point, because
    // condition_number() calls BackSolve() and BackSolve() returns immediately
    // when m_isFactored is false. The original ordering meant every reported
    // condition number was computed from an all-zero solution vector.)
    if (imp->m_print_cn)
    {
        double c = condition_number();
        feLog("\tcondition number (est.) ................... : %lg\n\n", c);
    }

    return true;
}

//-----------------------------------------------------------------------------
// Numeric factorization of the current (possibly scaled) matrix.
bool AccelerateSparseSolver::FactorNow()
{
    if (__builtin_available(macOS 10.13, *)) {

        if (imp->m_isFactored) {
            SparseCleanup(imp->ASF);
            imp->m_isFactored = false;
        }

        imp->ASF = SparseFactor(imp->ASS, imp->A, imp->SNFO);
        if (imp->ASF.status != SparseStatusOK) {
            fprintf(stderr, "\n\tERROR during factorization:");
            switch (imp->ASF.status) {
                case SparseFactorizationFailed:
                    fprintf(stderr, "\n\tFactorization failed!");
                    break;
                case SparseMatrixIsSingular:
                    fprintf(stderr, "\n\tSparse matrix is singular!");
                    break;
                case SparseInternalError:
                    fprintf(stderr, "\n\tSolver called internal error!");
                    break;
                case SparseParameterError:
                    fprintf(stderr, "\n\tSolver called parameter error!");
                    break;
                default:
                    fprintf(stderr, "\n\tUnknown error!");
            }
            return false;
        }

        imp->m_isFactored = true;
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------
bool AccelerateSparseSolver::BackSolve(double* x, double* b)
{
    // make sure we have work to do
    if (imp->m_pA->Rows() == 0) return true;

    const int strat = imp->strategy();
    if ((strat != SSIterative) && (imp->m_isFactored == false)) return true;

    const int  n      = imp->m_n;
    const bool bscale = (imp->m_scaling == PartitionScalingOn) && ((int)imp->m_scale.size() == n);

    // We solve A x = b by solving (D A D) y = D b and recovering x = D y.
    //
    // The previous version built the scaled right-hand side and then never used
    // it (B.data was still set to the unscaled b), and never undid the column
    // scaling on the solution. Both halves of the transformation have to be
    // applied or the answer is wrong by a factor of D on every entry.
    double* rhs = b;
    if (bscale) {
        imp->m_bs.resize(n);
        for (int i = 0; i < n; ++i) imp->m_bs[i] = imp->m_scale[i] * b[i];
        rhs = imp->m_bs.data();
    }

    DenseVector_Double X, B;
    X.count = B.count = n;
    X.data = x; B.data = rhs;

    // solve the system
    if (__builtin_available(macOS 10.13, *)) {
        if (strat == SSDirect) {
            SparseSolve(imp->ASF, B, X);
        }
        else if (strat == SSFactorPreconditioned) {
            // GMRES preconditioned by the stored (stale) factorization.
            SparseGMRESOptions gmres_opt = {};
            gmres_opt.maxIterations = imp->m_maxiter;
            gmres_opt.reportError   = MyReportError;
            gmres_opt.reportStatus  = (imp->m_print_level > 0) ? MyReportStatus : nullptr;
            gmres_opt.nvec          = imp->m_nvec;
            gmres_opt.rtol          = imp->m_rtol;
            gmres_opt.atol          = imp->m_atol;
            switch (imp->m_itrmthd) {
                case ITSparseDQGMRES: gmres_opt.variant = SparseVariantDQGMRES; break;
                case ITSparseGMRES:   gmres_opt.variant = SparseVariantGMRES;   break;
                // FGMRES is the safe default here: a factorization-based
                // preconditioner is fixed within a solve, but FGMRES also
                // tolerates it changing between solves as we refactor.
                default:              gmres_opt.variant = SparseVariantFGMRES;  break;
            }

            SparseOpaquePreconditioner_Double pc = {};
            pc.type  = SparsePreconditionerUser;
            pc.mem   = imp;
            pc.apply = ApplyStoredFactor;

            imp->SStatus = SparseSolve(SparseGMRES(gmres_opt), imp->A, B, X, pc);
            imp->m_pcAge++;

            // If the stale factorization was no longer good enough, refresh it
            // and retry once. This mirrors the nonlinear strategy in the paper:
            // when the quasi-Newton updates stop working, do another direct
            // solve rather than pressing on.
            if ((imp->SStatus != SparseIterativeConverged) && (imp->m_pcAge > 1)) {
                if (imp->m_print_level > 0)
                    feLog("\tGMRES failed with a preconditioner of age %d; refactoring and retrying\n", imp->m_pcAge);
                if (FactorNow()) {
                    imp->m_pcAge      = 1;
                    imp->m_pcRefactor = false;
                    imp->SStatus = SparseSolve(SparseGMRES(gmres_opt), imp->A, B, X, pc);
                }
            }

            // Whatever happened, make sure the next Factor() refreshes if we
            // still did not converge.
            if (imp->SStatus != SparseIterativeConverged) imp->m_pcRefactor = true;
        }
        else {
            // Pick the preconditioner once, rather than hard-wiring a different
            // one per GMRES variant as the previous version did.
            auto pc = SparsePreconditionerDiagonal;
            switch (imp->m_precond) {
                case PCNone:        pc = SparsePreconditionerNone;        break;
                case PCDiagonal:    pc = SparsePreconditionerDiagonal;    break;
                case PCDiagScaling: pc = SparsePreconditionerDiagScaling; break;
                default:            pc = SparsePreconditionerDiagonal;    break;
            }

            switch (imp->m_itrmthd) {
                case ITSparseConjugateGradient:
                {
                    // Braces are required: a case label may not jump over a
                    // variable with an initializer, and value-initializing the
                    // options struct is what guarantees the fields we do not
                    // set explicitly are zero rather than stack garbage.
                    SparseCGOptions cg_opt = {};
                    cg_opt.maxIterations = imp->m_maxiter;
                    cg_opt.rtol          = imp->m_rtol;
                    cg_opt.atol          = imp->m_atol;
                    cg_opt.reportError   = MyReportError;
                    cg_opt.reportStatus  = (imp->m_print_level > 0) ? MyReportStatus : nullptr;
                    imp->SStatus = SparseSolve(SparseConjugateGradient(cg_opt), imp->A, B, X, pc);
                    break;
                }
                case ITSparseGMRES:
                case ITSparseDQGMRES:
                case ITSparseFGMRES:
                {
                    SparseGMRESOptions gmres_opt = {};
                    gmres_opt.maxIterations = imp->m_maxiter;
                    gmres_opt.reportError   = MyReportError;
                    gmres_opt.reportStatus  = (imp->m_print_level > 0) ? MyReportStatus : nullptr;
                    gmres_opt.nvec          = imp->m_nvec;   // <= 0 lets Accelerate use its default of 16
                    gmres_opt.rtol          = imp->m_rtol;
                    gmres_opt.atol          = imp->m_atol;
                    switch (imp->m_itrmthd) {
                        case ITSparseGMRES:   gmres_opt.variant = SparseVariantGMRES;   break;
                        case ITSparseDQGMRES: gmres_opt.variant = SparseVariantDQGMRES; break;
                        case ITSparseFGMRES:  gmres_opt.variant = SparseVariantFGMRES;  break;
                    }
                    imp->SStatus = SparseSolve(SparseGMRES(gmres_opt), imp->A, B, X, pc);
                    break;
                }
                case ITSparseLSMR:
                {
                    SparseLSMROptions lsmr_opt = {};
                    lsmr_opt.rtol         = imp->m_rtol;
                    lsmr_opt.atol         = imp->m_atol;
                    lsmr_opt.reportError  = MyReportError;
                    lsmr_opt.reportStatus = (imp->m_print_level > 0) ? MyReportStatus : nullptr;
                    imp->SStatus = SparseSolve(SparseLSMR(lsmr_opt), imp->A, B, X, SparsePreconditionerDiagScaling);
                    break;
                }
                default:
                    fprintf(stderr, "\nUnknown iterative method (%d)!\n", imp->m_itrmthd);
                    return false;
            }
        }

        // Status reporting is shared by both iterative strategies.
        if ((strat != SSDirect) && (imp->SStatus != SparseIterativeConverged)) {
            fprintf(stderr, "\n\tERROR during iterative solution:\n");
            switch (imp->SStatus) {
                case SparseIterativeIllConditioned:
                    fprintf(stderr, "\n\tIll-conditioned system!\n");
                    break;
                case SparseIterativeInternalError:
                    fprintf(stderr, "\n\tInternal failure!\n");
                    break;
                case SparseIterativeMaxIterations:
                    fprintf(stderr, "\n\tExceeded maximum iteration limit!\n");
                    break;
                case SparseIterativeParameterError:
                    fprintf(stderr, "\n\tError with one or more parameters!\n");
                    break;
                default:
                    fprintf(stderr, "\n\tUnknown error!\n");
            }
            return false;
        }
    } else {
        // Fallback on earlier versions
        fprintf(stderr, "\nERROR during back solve: ");
        fprintf(stderr, "\nAccelerate solver not available for macOS earlier than 10.13!");
        return false;
    }

    // Undo the column scaling: x = D y.
    if (bscale) {
        for (int i = 0; i < n; ++i) x[i] *= imp->m_scale[i];
    }

    // update stats
    UpdateStats(1);

    return true;
}

//-----------------------------------------------------------------------------
// Partition-aware diagonal scaling.
//
// Builds a diagonal matrix D with one constant factor per FEBio partition and
// forms A' = D A D. For the velocity-dilatation CFD formulation the tangent has
// two partitions whose entries differ by roughly the fluid bulk modulus, so a
// single factor per partition removes essentially all of the imbalance.
//
// D is symmetric, so a symmetric A stays symmetric and the Cholesky/LDLT path
// remains valid.
//
// The scaled values go into a separate array; the assembled FEBio matrix is
// never modified, so the caller can still read the unscaled operator (as
// condition_number() does).
bool AccelerateSparseSolver::ComputeAndApplyScaling()
{
    const int n   = imp->m_n;
    const int nnz = imp->m_nnz;

    imp->m_scale.assign(n, 1.0);

    if (imp->m_scaling != PartitionScalingOn) {
        imp->m_scaledValues.clear();
        return true;
    }

    imp->m_scaledValues.resize(nnz);

    double* values   = imp->m_pA->Values();
    int*    pointers = imp->m_pA->Pointers();  // column starts, size ncols+1, 0-based
    int*    indices  = imp->m_pA->Indices();   // row indices per column

    // --- Determine partition boundaries from FEBio's own block structure ---
    const int nparts = Partitions();

    std::vector<int> partOf(n, 0);
    int neff = 1;

    if (nparts > 1) {
        std::vector<int> partStart(nparts + 1, 0);
        for (int p = 0; p < nparts; ++p)
            partStart[p + 1] = partStart[p] + GetPartitionSize(p);

        // The partitions must tile the equation set exactly. If they do not,
        // fall back to a single partition rather than indexing out of range --
        // the previous version indexed scalePart[] with a partition id that
        // could exceed its size whenever Partitions() returned 0.
        if (partStart[nparts] != n) {
            fprintf(stderr, "\nWARNING: partition sizes (%d) do not sum to the number of equations (%d);"
                            " partition scaling disabled.\n", partStart[nparts], n);
        }
        else {
            for (int p = 0; p < nparts; ++p)
                for (int i = partStart[p]; i < partStart[p + 1]; ++i)
                    partOf[i] = p;
            neff = nparts;
        }
    }

    // --- Largest |A_ij| touching each partition, by row and by column ---
    std::vector<double> maxAbsPart(neff, 0.0);

    const int ncols = imp->m_pA->Columns();
    for (int col = 0; col < ncols; ++col)
    {
        const int colPart = partOf[col];
        for (int idx = pointers[col]; idx < pointers[col + 1]; ++idx)
        {
            const int    row     = indices[idx];
            const double v       = std::abs(values[idx]);
            const int    rowPart = partOf[row];
            if (v > maxAbsPart[rowPart]) maxAbsPart[rowPart] = v;
            if (v > maxAbsPart[colPart]) maxAbsPart[colPart] = v;
        }
    }

    // --- One scale factor per partition: D_p = 1/sqrt(max|A|_p) ---
    std::vector<double> scalePart(neff, 1.0);
    for (int p = 0; p < neff; ++p)
        if (maxAbsPart[p] > 1e-30)
            scalePart[p] = 1.0 / std::sqrt(maxAbsPart[p]);

    for (int i = 0; i < n; ++i)
        imp->m_scale[i] = scalePart[partOf[i]];

    // --- Apply A'_ij = D_i * A_ij * D_j ---
    for (int col = 0; col < ncols; ++col)
    {
        const double sCol = imp->m_scale[col];
        for (int idx = pointers[col]; idx < pointers[col + 1]; ++idx)
        {
            const int row = indices[idx];
            imp->m_scaledValues[idx] = values[idx] * imp->m_scale[row] * sCol;
        }
    }

    if (imp->m_print_level > 0) {
        for (int p = 0; p < neff; ++p)
            feLog("\tpartition %d scaling ....................... : %lg\n", p, scalePart[p]);
    }

    return true;
}

//-----------------------------------------------------------------------------
// This algorithm (naively) estimates the condition number. It is based on the observation that
// for a linear system of equations A.x = b, the following holds
// || A^-1 || >= ||x||.||b||
// Thus the condition number can be estimated by
// c = ||A||.||A^-1|| >= ||A|| . ||x|| / ||b||
// This algorithm tries for some random b vectors with norm ||b||=1 to maxize the ||x||.
// The returned value will be an underestimate of the condition number.
//
// Note this reports the condition number of the UNSCALED operator, because
// BackSolve() undoes the scaling internally and infNorm() is taken on the
// assembled matrix. Comparing this against a run with partition_scaling on is
// therefore not meaningful; instrument ComputeAndApplyScaling() instead.
double AccelerateSparseSolver::condition_number()
{
    // This assumes that the factorization is already done!
    int N = imp->m_pA->Rows();
    if (N == 0) return 0.0;

    // get the norm of the matrix
    double normA = imp->m_pA->infNorm();

    // estimate the norm of the inverse of A
    double normAi = 0.0;

    // choose max iterations
    int iters = (N < 50 ? N : 50);

    std::vector<double> b(N, 0), x(N, 0);
    for (int i = 0; i < iters; ++i)
    {
        // create a random vector
        NumCore::randomVector(b, -1.0, 1.0);
        for (int j = 0; j < N; ++j) b[j] = (b[j] >= 0.0 ? 1.0 : -1.0);

        // calculate solution
        BackSolve(&x[0], &b[0]);

        double normx = NumCore::infNorm(x);
        if (normx > normAi) normAi = normx;

        // guard against division by zero when iters == 1
        int pct = (iters > 1) ? (100 * i) / (iters - 1) : 100;
        fprintf(stderr, "calculating condition number: %d%%\r", pct);
    }

    double c = normA*normAi;
    return c;
}

//-----------------------------------------------------------------------------
void AccelerateSparseSolver::Destroy()
{
    // Clean the symbolic and numeric factorizations independently. The previous
    // version cleaned neither unless m_isFactored was true, so a symbolic
    // factorization followed by a failed numeric factorization leaked.
    if (__builtin_available(macOS 10.13, *)) {
        if (imp->m_isFactored) SparseCleanup(imp->ASF);
        if (imp->m_isSymbolic) SparseCleanup(imp->ASS);
    }
    imp->m_isFactored = false;
    imp->m_isSymbolic = false;
    imp->m_pcAge      = 0;
    imp->m_pcRefactor = true;

    imp->m_scale.clear();
    imp->m_scaledValues.clear();
    imp->m_bs.clear();
}

#else
BEGIN_FECORE_CLASS(AccelerateSparseSolver, LinearSolver)
END_FECORE_CLASS();

AccelerateSparseSolver::AccelerateSparseSolver(FEModel* fem) : LinearSolver(fem), imp(nullptr) {}
AccelerateSparseSolver::~AccelerateSparseSolver() {}
bool AccelerateSparseSolver::PreProcess() { return false; }
bool AccelerateSparseSolver::Factor() { return false; }
bool AccelerateSparseSolver::BackSolve(double* x, double* y) { return false; }
void AccelerateSparseSolver::Destroy() {}
SparseMatrix* AccelerateSparseSolver::CreateSparseMatrix(Matrix_Type ntype) { return nullptr; }
bool AccelerateSparseSolver::SetSparseMatrix(SparseMatrix* pA) { return false; }
void AccelerateSparseSolver::PrintConditionNumber(bool b) {}
double AccelerateSparseSolver::condition_number() { return 0; }
void AccelerateSparseSolver::UseIterativeFactorization(bool b) {}
void AccelerateSparseSolver::SetFactorizationType(int ftype) {}
void AccelerateSparseSolver::SetOrderMethod(int order) {}
void AccelerateSparseSolver::SetIterativeMethod(int method) {}
void AccelerateSparseSolver::SetPartitionScaling(PartitionScaling s) {}
void AccelerateSparseSolver::SetSolveStrategy(int strategy) {}
void AccelerateSparseSolver::SetPrintLevel(int printlevel) {}
bool AccelerateSparseSolver::IsIterative() const { return false; }
bool AccelerateSparseSolver::ComputeAndApplyScaling() { return false; }
bool AccelerateSparseSolver::FactorNow() { return false; }
bool AccelerateSparseSolver::ResolveParameters() { return false; }
#endif
