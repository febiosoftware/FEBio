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
#include "PastixSparseSolver.h"
#include "MatrixTools.h"
#include <FECore/log.h>
#include <cstring>
#include <iostream>
#include <pastix.h>
#include <spm.h>
#include <cmath>

PastixSparseSolver::PastixSparseSolver(FEModel* fem)
    : LinearSolver(fem)
{
}

PastixSparseSolver::~PastixSparseSolver()
{
    Destroy();
}

SparseMatrix* PastixSparseSolver::CreateSparseMatrix(Matrix_Type ntype)
{
    switch (ntype)
    {
        case REAL_SYMMETRIC:    m_mtype = SpmSymmetric; m_pA = new CompactSymmMatrix(1); break;
        case REAL_UNSYMMETRIC:
        case REAL_SYMM_STRUCTURE: m_mtype = SpmGeneral; m_pA = new CRSSparseMatrix(1);
            break;
        default:
            m_pA = new CRSSparseMatrix(1);
            break;
    }

    return m_pA;
}

bool PastixSparseSolver::SetSparseMatrix(SparseMatrix* pA)
{
    if (m_pA && m_factored) Destroy();
    m_pA = dynamic_cast<CompactMatrix*>(pA);
    m_mtype = SpmSymmetric;
    if (dynamic_cast<CRSSparseMatrix*>(pA)) m_mtype = SpmGeneral;
    return (m_pA != nullptr);
}

bool PastixSparseSolver::PreProcess()
{
    if (m_pA == nullptr) return false;

    m_n = m_pA->Rows();
    m_nrhs = 1;

    pastixInitParam(m_iparm, m_dparm);

//    m_iparm[IPARM_MODIFY_PARAMETER] = 1;
    m_iparm[IPARM_FLOAT]            = PastixDouble;
    m_iparm[IPARM_MTX_TYPE]         = m_mtype;
    m_iparm[IPARM_FACTORIZATION]    = (m_mtype == SpmSymmetric) ? PastixFactLDLT : PastixFactLU;
    m_iparm[IPARM_ITERMAX]          = 100;
    
    // new stuff recommended by Claude:
    
    // --- Ordering ---
    // Choose one; benchmark both on your real meshes
    m_iparm[IPARM_ORDERING] = PastixOrderMetis;   // or PastixOrderScotch
    
    // --- Static pivoting / numerical stability ---
    m_iparm[IPARM_STATIC_PIVOTING] = 1;              // enable
    m_dparm[DPARM_EPSILON_MAGN_CTRL] = 1e-31;        // magnitude control threshold
    // (tune per problem; smaller = stricter)
    
    // --- Iterative refinement after the direct solve ---
    m_iparm[IPARM_ITERMAX] = 20;                      // max refinement iterations
    m_dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;        // convergence tolerance
    // iparm[IPARM_REFINEMENT] can select the refinement algorithm,
    // e.g. PastixRefineGMRES (default) vs PastixRefineBiCGSTAB, etc.
    m_iparm[IPARM_REFINEMENT] = PastixRefineGMRES;
    
    // --- Threading / scheduler ---
    m_iparm[IPARM_THREAD_NBR] = 18;                    // set explicitly, don't rely on default
    m_iparm[IPARM_SCHEDULER]  = PastixSchedDynamic;   // or PastixSchedStatic, PastixSchedSequential,
    // PastixSchedParsec, PastixSchedStarPU
    // (Parsec/StarPU require the runtime to be built in)
    
    // --- Block low-rank compression (memory/speed for large 3D problems) ---
    m_iparm[IPARM_COMPRESS_WHEN]       = PastixCompressWhenDuring; // or PastixCompressBegin / PastixCompressNever
    m_iparm[IPARM_COMPRESS_METHOD]     = PastixCompressMethodPQRCP; // or PastixCompressMethodSVD
    m_iparm[IPARM_COMPRESS_MIN_WIDTH]  = 128;   // min supernode width to consider compressing
    m_iparm[IPARM_COMPRESS_MIN_HEIGHT] = 20;    // min off-diagonal block height to compress
    m_dparm[DPARM_COMPRESS_TOLERANCE]  = 1e-8;  // compression accuracy tradeoff
    
    // --- Verbosity (useful while tuning) ---
    m_iparm[IPARM_VERBOSE] = PastixVerboseYes;  // or PastixVerboseNo / PastixVerboseNot
    
    // Then proceed as usual:
    pastix_data_t *pastix_data = nullptr;
    pastixInit(&m_pastix, MPI_COMM_WORLD, m_iparm, m_dparm);

    m_initialized = true;
    m_factored = false;

    return LinearSolver::PreProcess();
}

bool PastixSparseSolver::CopyPatternToPastix()
{
    if (m_pA == nullptr) return false;

    const int n = m_pA->Rows();
    const int nnz = static_cast<int>(m_pA->NonZeroes());

    int* p = m_pA->Pointers();
    int* i = m_pA->Indices();

    if (p == nullptr || i == nullptr) return false;

    m_colptr.resize(n + 1);
    m_rowind.resize(nnz);
    m_perm.assign(n, 0);
    m_invp.assign(n, 0);

    const int offset = m_pA->Offset();

    for (int k = 0; k <= n; ++k)
        m_colptr[k] = static_cast<pastix_int_t>(p[k] - offset);

    for (int k = 0; k < nnz; ++k)
        m_rowind[k] = static_cast<pastix_int_t>(i[k] - offset);

    return true;
}

bool PastixSparseSolver::Factor()
{
    if (m_pA == nullptr) return false;
    if (m_pA->Rows() == 0) return true;
    
    if (!m_initialized) {
        if (!PreProcess()) return false;
    }
    
    if (!CopyPatternToPastix()) return false;
    
    // NEW: compute diagonal scaling and produce a scaled copy of the values
    // (matrix values change every Newton reformation, so recompute each Factor() call)
    if (!ComputeAndApplyScaling()) return false;
    
    std::vector<double> dummy_rhs(m_n, 0.0);
    
    m_iparm[IPARM_START_TASK] = PastixTaskOrdering;
    m_iparm[IPARM_END_TASK]   = PastixTaskNumfact;
    
    int ierr = pastix(
                      &m_pastix,
                      0,
                      static_cast<pastix_int_t>(m_n),
                      m_colptr.data(),
                      m_rowind.data(),
                      m_scaledValues.data(),      // CHANGED: was m_pA->Values()
                      m_perm.data(),
                      m_invp.data(),
                      dummy_rhs.data(),
                      static_cast<pastix_int_t>(m_nrhs),
                      m_iparm,
                      m_dparm
                      );
    
    if (ierr != 0) {
        std::cerr << "PaStiX factorization error: " << ierr << std::endl;
        return false;
    }
    
    m_factored = true;
    return true;
}
bool PastixSparseSolver::BackSolve(double* x, double* b)
{
    if (m_pA == nullptr) return false;
    if (m_pA->Rows() == 0) return true;
    if (!m_factored) return false;
    
    std::vector<double> rhs(m_n);
    for (int i = 0; i < m_n; ++i)
        rhs[i] = m_scale[i] * b[i];
    
    m_iparm[IPARM_START_TASK] = PastixTaskSolve;
    m_iparm[IPARM_END_TASK]   = PastixTaskRefine;
    m_iparm[IPARM_TRANSPOSE_SOLVE] = PastixTrans;   // NEW: compensate for CRSSparseMatrix's
    // row-major storage being fed to PaStiX's
    // column-major colptr/rowind parameters
    
    int ierr = pastix(
                      &m_pastix,
                      0,
                      static_cast<pastix_int_t>(m_n),
                      m_colptr.data(),
                      m_rowind.data(),
                      m_scaledValues.data(),
                      m_perm.data(),
                      m_invp.data(),
                      rhs.data(),
                      static_cast<pastix_int_t>(m_nrhs),
                      m_iparm,
                      m_dparm
                      );
    
    if (ierr != 0) {
        std::cerr << "PaStiX solve error: " << ierr << std::endl;
        return false;
    }
    
    for (int i = 0; i < m_n; ++i)
        x[i] = m_scale[i] * rhs[i];
    
    UpdateStats(1);
    return true;
}

void PastixSparseSolver::Destroy()
{
    if (m_pastix != nullptr) {
        pastixFinalize(&m_pastix);
    }

    m_initialized = false;
    m_factored = false;

    m_colptr.clear();
    m_rowind.clear();
    m_perm.clear();
    m_invp.clear();
}

// Computes symmetric diagonal (Jacobi-type) scaling: scale[i] = 1/sqrt(|A_ii|),
// and produces a scaled copy of the matrix values: A'_ij = scale[i] * A_ij * scale[j].
// This compensates for PaStiX 6.4's lack of an MC64-equivalent matching/scaling step.
bool PastixSparseSolver::ComputeAndApplyScaling()
{
    const int n   = m_n;
    const int nnz = static_cast<int>(m_rowind.size());
    
    m_scale.assign(n, 1.0);
    m_scaledValues.resize(nnz);
    
    double* values = m_pA->Values();
    
    // --- Step 1: determine partition boundaries from FEBio's own block structure ---
    // (m_part is inherited from LinearSolver; set via SetPartitions() before Factor())
    int nparts = Partitions();
    std::vector<int> partStart(nparts + 1, 0);
    for (int p = 0; p < nparts; ++p)
        partStart[p + 1] = partStart[p] + GetPartitionSize(p);
    
    // Map each equation index to its partition id
    std::vector<int> partOf(n, 0);
    for (int p = 0; p < nparts; ++p)
        for (int i = partStart[p]; i < partStart[p + 1]; ++i)
            partOf[i] = p;
    
    // --- Step 2: find max |A_ij| touching each row/column, accumulated per partition ---
    std::vector<double> maxAbsPart(nparts, 0.0);
    
    for (int col = 0; col < n; ++col)
    {
        int colPart = partOf[col];
        for (pastix_int_t idx = m_colptr[col]; idx < m_colptr[col + 1]; ++idx)
        {
            int row = static_cast<int>(m_rowind[idx]);
            double v = std::abs(values[idx]);
            int rowPart = partOf[row];
            if (v > maxAbsPart[rowPart]) maxAbsPart[rowPart] = v;
            if (v > maxAbsPart[colPart]) maxAbsPart[colPart] = v;
        }
    }
    
    // --- Step 3: one scale factor per partition (not per equation) ---
    std::vector<double> scalePart(nparts, 1.0);
    for (int p = 0; p < nparts; ++p)
    {
        if (maxAbsPart[p] > 1e-30)
            scalePart[p] = 1.0 / std::sqrt(maxAbsPart[p]);
    }
    
    for (int i = 0; i < n; ++i)
        m_scale[i] = scalePart[partOf[i]];
    
    // --- Step 4: apply A'_ij = scale[i] * A_ij * scale[j] ---
    for (int col = 0; col < n; ++col)
    {
        double sCol = m_scale[col];
        for (pastix_int_t idx = m_colptr[col]; idx < m_colptr[col + 1]; ++idx)
        {
            int row = static_cast<int>(m_rowind[idx]);
            m_scaledValues[idx] = values[idx] * m_scale[row] * sCol;
        }
    }
    
    return true;
}
