#include "stdafx.h"
#include "FEExplicitFluidSolver.h"
#include "FEFluidDomain3D.h"
#include "FEFluid.h"
#include "FEFluidAnalysis.h"

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FENode.h>
#include <FECore/FEAnalysis.h>
#include <FECore/log.h>

#include <algorithm>
#include <cmath>
#include <limits>


BEGIN_FECORE_CLASS(FEExplicitFluidSolver, FEFluidSolver)
    ADD_PARAMETER(m_dynDamping,    "dyn_damping");
    ADD_PARAMETER(m_capacityFloor, "capacity_floor");
END_FECORE_CLASS();

FEExplicitFluidSolver::FEExplicitFluidSolver(FEModel* fem)
    : FEFluidSolver(fem)
{
    m_dynDamping = 1.0;
    m_capacityFloor = 1.0e-30;

    // In this FEBio branch, rhoi=-1 selects alpha_f=alpha_m=gamma=1.
    // UpdateKinematics then reconstructs accepted rates from first differences.
    m_rhoi = -1.0;
    m_pred = 0;

    // The explicit path does not reform or factor a stiffness matrix.
    m_persistMatrix = true;
}

bool FEExplicitFluidSolver::Init()
{
    if (FEFluidSolver::Init() == false) return false;

    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    if (fem.GetCurrentStep()->m_nanalysis != FEFluidAnalysis::DYNAMIC)
    {
        feLogError("The RK4 explicit fluid solver requires a dynamic fluid analysis.");
        return false;
    }

    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEFluidDomain3D* dom = dynamic_cast<FEFluidDomain3D*>(&mesh.Domain(i));
        if (dom == nullptr)
        {
            feLogError("The RK4 explicit fluid solver currently supports FEFluidDomain3D domains only.");
            return false;
        }
        dom->SetTransientAnalysis();
    }

    m_capacity.assign(m_neq, 0.0);


    return true;
}

void FEExplicitFluidSolver::ZeroPrimaryRates()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        node.set(m_dofAW[0], 0.0);
        node.set(m_dofAW[1], 0.0);
        node.set(m_dofAW[2], 0.0);
        node.set(m_dofAEF, 0.0);
    }
}

bool FEExplicitFluidSolver::BuildLumpedOperators()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    m_capacity.assign(m_neq, 0.0);

    for (int nd = 0; nd < mesh.Domains(); ++nd)
    {
        FEFluidDomain3D* dom = dynamic_cast<FEFluidDomain3D*>(&mesh.Domain(nd));
        if (dom == nullptr) return false;

        FEFluid* mat = dynamic_cast<FEFluid*>(dom->GetMaterial());
        if (mat == nullptr)
        {
            feLogError("A fluid domain does not reference an FEFluid material.");
            return false;
        }

        const int NE = dom->Elements();
        int failed = 0;

#pragma omp parallel for schedule(static) reduction(|:failed)
        for (int iel = 0; iel < NE; ++iel)
        {
            FESolidElement& el = dom->Element(iel);
            const int neln = el.Nodes();

            vector<double> mv(neln, 0.0);
            vector<double> ce(neln, 0.0);

            bool elementOK = true;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEFluidMaterialPoint* fp = mp.ExtractData<FEFluidMaterialPoint>();
                if (fp == nullptr)
                {
                    elementOK = false;
                    break;
                }

                const double J = 1.0 + fp->m_ef;
                const double rho = mat->Density(mp);
                const double dv = dom->detJ0(el, n) * el.GaussWeights()[n];
                const double* H = el.H(n);

                if (!(J > 0.0) || !std::isfinite(J) ||
                    !(rho > 0.0) || !std::isfinite(rho) ||
                    !(dv > 0.0) || !std::isfinite(dv))
                {
                    elementOK = false;
                    break;
                }

                for (int a = 0; a < neln; ++a)
                {
                    mv[a] += rho * H[a] * dv;
                    ce[a] += H[a] * dv / J;
                }
            }

            if (!elementOK)
            {
                failed = 1;
                continue;
            }

            for (int a = 0; a < neln; ++a)
            {
                FENode& node = mesh.Node(el.m_node[a]);

                for (int k = 0; k < 3; ++k)
                {
                    const int eq = node.m_ID[m_dofW[k]];
                    if ((eq >= 0) && (eq < m_neq))
                    {
#pragma omp atomic
                        m_capacity[eq] += mv[a];
                    }
                }

                const int eqe = node.m_ID[m_dofEF[0]];
                if ((eqe >= 0) && (eqe < m_neq))
                {
#pragma omp atomic
                    m_capacity[eqe] += ce[a];
                }
            }
        }

        if (failed)
        {
            feLogError("Invalid material-point data while assembling RK4 capacities in domain %d.", nd + 1);
            return false;
        }
    }

    int active = 0;
    double cmin = std::numeric_limits<double>::max();
    double cmax = 0.0;

#pragma omp parallel for reduction(+:active) reduction(min:cmin) reduction(max:cmax)
    for (int i = 0; i < m_neq; ++i)
    {
        const double c = m_capacity[i];
        if (c > m_capacityFloor)
        {
            cmin = std::min(cmin, c);
            cmax = std::max(cmax, c);
            active += 1;
        }
    }

    if (active == 0)
    {
        feLogError("No positive RK4 fluid capacities were assembled.");
        return false;
    }

    return true;
}

void FEExplicitFluidSolver::InsertPrescribedIncrements(
    const vector<double>& prescribed, vector<double>& increment, double scale) const
{
    const FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        const FENode& node = mesh.Node(i);

        for (int k = 0; k < 3; ++k)
        {
            const int id = node.m_ID[m_dofW[k]];
            if (id < -1)
            {
                const int eq = -id - 2;
                if ((eq >= 0) && (eq < (int)increment.size()))
                    increment[eq] = scale * prescribed[eq];
            }
        }

        const int id = node.m_ID[m_dofEF[0]];
        if (id < -1)
        {
            const int eq = -id - 2;
            if ((eq >= 0) && (eq < (int)increment.size()))
                increment[eq] = scale * prescribed[eq];
        }
    }
}

bool FEExplicitFluidSolver::FormExplicitRate(
    const vector<double>& R, vector<double>& rate) const
{
    rate.assign(m_neq, 0.0);

    int active = 0;
    int invalidEq = m_neq;

#pragma omp parallel for reduction(+:active) reduction(min:invalidEq)
    for (int i = 0; i < m_neq; ++i)
    {
        if (!std::isfinite(R[i]))
        {
            invalidEq = std::min(invalidEq, i);
            continue;
        }

        const double c = m_capacity[i];
        if (c > m_capacityFloor)
        {
            const double ri = m_dynDamping * R[i] / c;
            if (!std::isfinite(ri))
            {
                invalidEq = std::min(invalidEq, i);
                continue;
            }
            rate[i] = ri;
            active += 1;
        }
    }

    if (invalidEq < m_neq)
    {
        feLogError("Non-finite residual or RK4 rate in fluid equation %d.", invalidEq);
        return false;
    }

    return (active > 0);
}

bool FEExplicitFluidSolver::EvaluateStage(
    const vector<double>& stageIncrement, vector<double>& rate)
{
    vector<double> u(stageIncrement);
    Update(u);

    ZeroPrimaryRates();
    UpdateModel();

    if (BuildLumpedOperators() == false) return false;

    vector<double> R(m_neq, 0.0);
    if (Residual(R) == false) return false;

    return FormExplicitRate(R, rate);
}

void FEExplicitFluidSolver::LogPrimaryState(const char* label) const
{
    const FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
    double vmax = 0.0;
    double emax = 0.0;
    double Jmin = std::numeric_limits<double>::max();

#pragma omp parallel for reduction(max:vmax,emax) reduction(min:Jmin)
    for (int i = 0; i < NN; ++i)
    {
        const FENode& node = mesh.Node(i);
        const vec3d v = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2]);
        const double e = node.get(m_dofEF[0]);
        const double J = 1.0 + e;

        vmax = std::max(vmax, v.norm());
        emax = std::max(emax, std::fabs(e));
        Jmin = std::min(Jmin, J);
    }

    feLog("%s: vmax=%.16g, emax=%.16g, minJ=%.16g\n",
        label, vmax, emax, Jmin);
}

bool FEExplicitFluidSolver::CheckExplicitState() const
{
    const FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
    int invalidNode = NN;

#pragma omp parallel for reduction(min:invalidNode)
    for (int i = 0; i < NN; ++i)
    {
        const FENode& node = mesh.Node(i);
        const double e = node.get(m_dofEF[0]);
        const double J = 1.0 + e;

        bool bad = (!std::isfinite(e) || !(J > 0.0));
        if (!bad)
        {
            for (int k = 0; k < 3; ++k)
            {
                if (!std::isfinite(node.get(m_dofW[k])))
                {
                    bad = true;
                    break;
                }
            }
        }

        if (bad) invalidNode = std::min(invalidNode, i);
    }

    if (invalidNode < NN)
    {
        const FENode& node = mesh.Node(invalidNode);
        const double J = 1.0 + node.get(m_dofEF[0]);
        if (!std::isfinite(J) || !(J > 0.0))
            feLogError("RK4 fluid update produced invalid J at node %d: J = %.16g", invalidNode + 1, J);
        else
            feLogError("RK4 fluid update produced non-finite velocity at node %d.", invalidNode + 1);
        return false;
    }

    return true;
}

bool FEExplicitFluidSolver::Quasin()
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    const double dt = tp.timeIncrement;

    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        feLogError("The RK4 explicit fluid solver requires a positive finite time increment.");
        return false;
    }

    m_niter = 0;
    m_nrhs = 0;
    m_nref = 0;
    m_ntotref = 0;

    PrepStep();
    const vector<double> prescribed = m_ui;

    vector<double> y0(m_neq, 0.0);
    vector<double> stage(m_neq, 0.0);
    vector<double> k1, k2, k3, k4;

    if (EvaluateStage(y0, k1) == false) return false;

#pragma omp parallel for
    for (int i = 0; i < m_neq; ++i) stage[i] = 0.5 * dt * k1[i];
    InsertPrescribedIncrements(prescribed, stage, 0.5);
    if (EvaluateStage(stage, k2) == false) return false;

#pragma omp parallel for
    for (int i = 0; i < m_neq; ++i) stage[i] = 0.5 * dt * k2[i];
    InsertPrescribedIncrements(prescribed, stage, 0.5);
    if (EvaluateStage(stage, k3) == false) return false;

#pragma omp parallel for
    for (int i = 0; i < m_neq; ++i) stage[i] = dt * k3[i];
    InsertPrescribedIncrements(prescribed, stage, 1.0);
    if (EvaluateStage(stage, k4) == false) return false;

    vector<double> du(m_neq, 0.0);
    double maxRate = 0.0;
    double maxDu = 0.0;

#pragma omp parallel for reduction(max:maxRate,maxDu)
    for (int i = 0; i < m_neq; ++i)
    {
        const double r = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
        du[i] = dt * r;
        maxRate = std::max(maxRate, std::fabs(r));
        maxDu = std::max(maxDu, std::fabs(du[i]));
    }

    InsertPrescribedIncrements(prescribed, du, 1.0);

    feLog("RK4 fluid update: max|rate|=%.16g, max|du|=%.16g\n", maxRate, maxDu);

    Update(du);
    LogPrimaryState("RK4 post-update state");

    UpdateIncrements(m_Ui, du, false);

    if (CheckExplicitState() == false)
    {
        feLogError("RK4 explicit state check failed.");
        return false;
    }

    UpdateIncrements(m_Ut, m_Ui, true);
    zero(m_Ui);

    m_niter = 1;
    fem.DoCallback(CB_MINOR_ITERS);
    return true;
}

void FEExplicitFluidSolver::Serialize(DumpStream& ar)
{
    FEFluidSolver::Serialize(ar);
    ar & m_capacity;
    ar & m_dynDamping & m_capacityFloor;
}
