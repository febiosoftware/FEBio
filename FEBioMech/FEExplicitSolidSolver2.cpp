#include "stdafx.h"
#include "FEExplicitSolidSolver2.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticShellDomain.h"
#include "FEElasticEASShellDomain.h"
#include "FEElasticANSShellDomain.h"
#include "FEElasticBeamDomain.h"
#include "FEElasticDomain.h"
#include "FEElasticMaterial.h"
#include "FESolidAnalysis.h"
#include "FEContactInterface.h"

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FENode.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FESurfacePairConstraint.h>
#include <FECore/FEModelLoad.h>
#include <FECore/FELinearSystem.h>
#include <FECore/log.h>/Users/gerard/GitHub/FEBio-safe/FEBioMech/FEBioMech.cpp

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
    inline bool IsFiniteVec3d(const vec3d& v)
    {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    }

    inline int EquationIndex(int id)
    {
        return (id >= 0 ? id : (id < -1 ? -id - 2 : -1));
    }

    // Intercepts element mass-matrix assembly and accumulates the row sums
    // directly, avoiding dependence on the sparse matrix implementation.
    class FERowSumMassSystem : public FELinearSystem
    {
    public:
        FERowSumMassSystem(FEModel* fem, FEGlobalMatrix& K,
            vector<double>& dummyF, vector<double>& dummyU,
            vector<double>& lumpedMass)
            : FELinearSystem(fem, K, dummyF, dummyU, false),
              m_lumpedMass(lumpedMass)
        {
        }

        void Assemble(const FEElementMatrix& ke) override
        {
            const int nr = ke.rows();
            const int nc = ke.columns();
            if ((nr == 0) || (nc == 0)) return;

            const vector<int>& lmi = ke.RowIndices();
            for (int i = 0; i < nr; ++i)
            {
                const int I = EquationIndex(lmi[i]);
                if ((I < 0) || (I >= (int)m_lumpedMass.size())) continue;

                double rowSum = 0.0;
                for (int j = 0; j < nc; ++j) rowSum += ke[i][j];

                if (std::isfinite(rowSum))
                {
#pragma omp atomic
                    m_lumpedMass[I] += rowSum;
                }
            }
        }

    private:
        vector<double>& m_lumpedMass;
    };
}

BEGIN_FECORE_CLASS(FEExplicitSolidSolver2, FESolidSolver2)
    ADD_PARAMETER(m_dynDamping, "dyn_damping");
    ADD_PARAMETER(m_massFloor,  "mass_floor");
    ADD_PARAMETER(m_minJ,       "min_jacobian");
END_FECORE_CLASS();

FEExplicitSolidSolver2::FEExplicitSolidSolver2(FEModel* fem)
    : FESolidSolver2(fem)
{
    m_dynDamping = 1.0;
    m_massFloor = 1.0e-30;
    m_minJ = 1.0e-12;
    m_init_accelerations = false;
    m_persistMatrix = true;
}

bool FEExplicitSolidSolver2::Init()
{
    if (FESolidSolver2::Init() == false) return false;

    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    if (fem.GetCurrentStep()->m_nanalysis != FESolidAnalysis::DYNAMIC)
    {
        feLogError("The explicit solid solver requires a dynamic analysis.");
        return false;
    }

    if (fem.NonlinearConstraints() != 0)
    {
        feLogError("The explicit solid solver does not support nonlinear-constraint equations.");
        return false;
    }

    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        if (mesh.Node(i).m_rid >= 0)
        {
            feLogError("The explicit solid solver does not support rigid bodies or rigid nodes.");
            return false;
        }
    }

    // Contact is supported in penalty form. Augmented-Lagrangian contact needs
    // an outer augmentation loop and is therefore rejected here.
    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        FEContactInterface* ci = dynamic_cast<FEContactInterface*>(spc);
        if (spc && spc->IsActive() && (ci == nullptr))
        {
            feLogError("Surface-pair constraint %d is not an FEContactInterface.", i + 1);
            return false;
        }
        if (ci && ci->IsActive() && (ci->m_laugon == FECore::AUGLAG_METHOD))
        {
            feLogError("The explicit solid solver currently supports penalty contact only; contact %d uses augmented Lagrangian enforcement.", i + 1);
            return false;
        }
    }

    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& d = mesh.Domain(i);
        FEElasticSolidDomain* solid = dynamic_cast<FEElasticSolidDomain*>(&d);
        FEElasticShellDomain* shell = dynamic_cast<FEElasticShellDomain*>(&d);
        FEElasticBeamDomain* beam = dynamic_cast<FEElasticBeamDomain*>(&d);

        if ((solid == nullptr) && (shell == nullptr) && (beam == nullptr))
        {
            feLogError("Domain %d is unsupported. Supported domain families are elastic solids, elastic shells, and elastic beams.", i + 1);
            return false;
        }

        if (solid) solid->SetDynamicUpdateFlag(true);
        if (shell) shell->SetDynamicUpdateFlag(true);
    }

    m_lumpedMass.assign(m_neq, 0.0);
    m_velocity.assign(m_neq, 0.0);
    m_acceleration.assign(m_neq, 0.0);
    m_stageVelocity.assign(m_neq, 0.0);
    m_stageAcceleration.assign(m_neq, 0.0);

    CaptureAcceptedVelocity();
    if (BuildLumpedMass() == false) return false;

    return true;
}

void FEExplicitSolidSolver2::PrepStep()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    FETimeInfo& tp = fem.GetTime();

    tp.augmentation = 0;
    tp.currentIteration = 0;
    zero(m_Ui);

    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_rp = node.m_rt;
        node.m_vp = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
        node.m_ap = node.m_at;
        node.m_dp = node.m_dt;
        node.UpdateValues();
    }

    zero(m_ui);
    for (int i = 0; i < fem.BoundaryConditions(); ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.PrepStep(m_ui);
    }

    fem.GetLinearConstraintManager().PrepStep();

    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEDomain& dom = mesh.Domain(i);
        if (dom.IsActive()) dom.PreSolveUpdate(tp);
    }

    UpdateModel();

    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        if (spc && spc->IsActive()) spc->PrepStep();
    }

    for (int i = 0; i < fem.ModelLoads(); ++i)
    {
        FEModelLoad* load = fem.ModelLoad(i);
        if (load && load->IsActive()) load->PrepStep();
    }

    m_baugment = false;
}

bool FEExplicitSolidSolver2::BuildLumpedMass()
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    if (m_pK == nullptr)
    {
        feLogError("The global matrix wrapper was not allocated before explicit mass assembly.");
        return false;
    }

    m_lumpedMass.assign(m_neq, 0.0);
    vector<double> dummyF(m_neq, 0.0), dummyU(m_neq, 0.0);
    FERowSumMassSystem LS(&fem, *m_pK, dummyF, dummyU, m_lumpedMass);

    // FEElasticDomain::MassMatrix dispatches to the domain-specific solid,
    // shell, or beam element mass implementation. The custom linear system
    // intercepts each element matrix and accumulates its row sums.
    for (int nd = 0; nd < mesh.Domains(); ++nd)
    {
        FEDomain& baseDomain = mesh.Domain(nd);
        FEElasticDomain* dom = dynamic_cast<FEElasticDomain*>(&baseDomain);
        if (dom == nullptr)
        {
            feLogError("Domain %d is not an FEElasticDomain.", nd + 1);
            return false;
        }
        if (baseDomain.IsActive()) dom->MassMatrix(LS, 1.0);
    }

    int active = 0;
    double mmin = std::numeric_limits<double>::max();
    double mmax = 0.0;
    int invalid = m_neq;
#pragma omp parallel for reduction(+:active) reduction(min:mmin,invalid) reduction(max:mmax)
    for (int i = 0; i < m_neq; ++i)
    {
        const double m = m_lumpedMass[i];
        if (!std::isfinite(m) || (m < -m_massFloor))
        {
            invalid = std::min(invalid, i);
        }
        else if (m > m_massFloor)
        {
            active += 1;
            mmin = std::min(mmin, m);
            mmax = std::max(mmax, m);
        }
    }

    if (invalid < m_neq)
    {
        feLogError("Invalid row-sum lumped mass in equation %d.", invalid);
        return false;
    }
    if (active == 0)
    {
        feLogError("No positive lumped masses were assembled.");
        return false;
    }

    feLog("Explicit solid/shell/beam lumped mass: active=%d, min=%.16g, max=%.16g\n", active, mmin, mmax);
    return true;
}

void FEExplicitSolidSolver2::GatherNodeRate(const FENode& node,
    const FEDofList& primary, const FEDofList& rate,
    vector<double>& globalRate) const
{
    const vec3d v = node.get_vec3d(rate[0], rate[1], rate[2]);
    const double a[3] = { v.x, v.y, v.z };
    for (int k = 0; k < 3; ++k)
    {
        const int eq = EquationIndex(node.m_ID[primary[k]]);
        if ((eq >= 0) && (eq < (int)globalRate.size())) globalRate[eq] = a[k];
    }
}

void FEExplicitSolidSolver2::CaptureAcceptedVelocity()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    std::fill(m_velocity.begin(), m_velocity.end(), 0.0);

    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        const FENode& node = mesh.Node(i);
        GatherNodeRate(node, m_dofU,  m_dofV,  m_velocity);
        GatherNodeRate(node, m_dofSU, m_dofSV, m_velocity);
        GatherNodeRate(node, m_dofQ,  m_dofBW, m_velocity);
    }
}

void FEExplicitSolidSolver2::InsertPrescribedPrimaryDofs(
    const vector<double>& prescribed, vector<double>& values,
    double velocityScale) const
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    const FEDofList* primary[3] = { &m_dofU, &m_dofSU, &m_dofQ };

    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        for (int family = 0; family < 3; ++family)
        {
            const FEDofList& dofs = *primary[family];
            for (int k = 0; k < 3; ++k)
            {
                const int id = node.m_ID[dofs[k]];
                if (id < -1)
                {
                    const int eq = -id - 2;
                    if ((eq >= 0) && (eq < (int)values.size()))
                        values[eq] = prescribed[eq] * velocityScale;
                }
            }
        }
    }
}

void FEExplicitSolidSolver2::SetStageNodeRate(FENode& node,
    const FEDofList& primary, const FEDofList& rate,
    const FEDofList* accelerationRate) const
{
    vec3d v(0, 0, 0), a(0, 0, 0);
    for (int k = 0; k < 3; ++k)
    {
        const int eq = EquationIndex(node.m_ID[primary[k]]);
        if ((eq >= 0) && (eq < (int)m_stageVelocity.size()))
        {
            const double vk = m_stageVelocity[eq];
            const double ak = m_stageAcceleration[eq];
            if (k == 0) { v.x = vk; a.x = ak; }
            else if (k == 1) { v.y = vk; a.y = ak; }
            else { v.z = vk; a.z = ak; }
        }
    }
    node.set_vec3d(rate[0], rate[1], rate[2], v);
    if (accelerationRate)
        node.set_vec3d((*accelerationRate)[0], (*accelerationRate)[1], (*accelerationRate)[2], a);
    else
        node.m_at = a;
}

void FEExplicitSolidSolver2::UpdateKinematics(vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    vector<double> U(m_Ut.size(), 0.0);
#pragma omp parallel for
    for (int i = 0; i < (int)U.size(); ++i) U[i] = ui[i] + m_Ui[i] + m_Ut[i];

    scatter3(U, mesh, m_dofU[0],  m_dofU[1],  m_dofU[2]);
    scatter3(U, mesh, m_dofSU[0], m_dofSU[1], m_dofSU[2]);
    scatter3(U, mesh, m_dofQ[0],  m_dofQ[1],  m_dofQ[2]);

    for (int i = 0; i < fem.BoundaryConditions(); ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.Update();
    }

    FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
    if (LCM.LinearConstraints() > 0) LCM.Update();

#pragma omp parallel for
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
        node.m_dt = node.m_d0
            + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2])
            - node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);

        SetStageNodeRate(node, m_dofU,  m_dofV,  nullptr);
        SetStageNodeRate(node, m_dofSU, m_dofSV, &m_dofSA);
        SetStageNodeRate(node, m_dofQ,  m_dofBW, &m_dofBA);
    }

    for (int i = 0; i < fem.SurfacePairConstraints(); ++i)
    {
        FESurfacePairConstraint* spc = fem.SurfacePairConstraint(i);
        if (spc && spc->IsActive()) spc->Update(m_Ui, ui);
    }
}

bool FEExplicitSolidSolver2::FormAcceleration(
    const vector<double>& R, vector<double>& acceleration) const
{
    acceleration.assign(m_neq, 0.0);
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

        const double mass = m_lumpedMass[i];
        if (mass > m_massFloor)
        {
            const double a = m_dynDamping * R[i] / mass;
            if (!std::isfinite(a)) invalidEq = std::min(invalidEq, i);
            else { acceleration[i] = a; active += 1; }
        }
    }

    if (invalidEq < m_neq)
    {
        feLogError("Non-finite residual or acceleration in equation %d.", invalidEq);
        return false;
    }
    if (active == 0)
    {
        feLogError("No active massive equations were found during the explicit update.");
        return false;
    }
    return true;
}

bool FEExplicitSolidSolver2::CheckElementJacobians(const char* stage, double& minJ) const
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    minJ = std::numeric_limits<double>::max();
    int badDomain = -1, badElement = -1, badPoint = -1;

    for (int nd = 0; nd < mesh.Domains(); ++nd)
    {
        FEDomain& baseDomain = mesh.Domain(nd);
        FEElasticDomain* dom = dynamic_cast<FEElasticDomain*>(&baseDomain);
        if (dom == nullptr) continue;

        for (int iel = 0; iel < baseDomain.Elements(); ++iel)
        {
            FEElement& el = baseDomain.ElementRef(iel);
            if (!el.isActive()) continue;
            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                FEElasticMaterialPoint* ep = mp.ExtractData<FEElasticMaterialPoint>();
                if (ep == nullptr) continue;
                const double J = ep->m_J;
                minJ = std::min(minJ, J);
                if (!std::isfinite(J) || !(J > m_minJ))
                {
                    badDomain = nd; badElement = iel; badPoint = n;
                    break;
                }
            }
            if (badDomain >= 0) break;
        }
        if (badDomain >= 0) break;
    }

    if (badDomain >= 0)
    {
        feLogError("Explicit solid %s state has invalid J in domain %d, element %d, integration point %d.",
            stage, badDomain + 1, badElement + 1, badPoint + 1);
        return false;
    }
    if (minJ == std::numeric_limits<double>::max()) minJ = 1.0;
    return true;
}

bool FEExplicitSolidSolver2::EvaluateAcceleration(
    const vector<double>& stageU, const vector<double>& stageV,
    vector<double>& acceleration)
{
    m_stageVelocity = stageV;
    m_stageAcceleration.assign(m_neq, 0.0);

    vector<double> u(stageU);
    Update(u);

    double minJ = 0.0;
    if (!CheckElementJacobians("trial", minJ)) return false;

    vector<double> R(m_neq, 0.0);
    if (Residual(R) == false) return false;
    return FormAcceleration(R, acceleration);
}

bool FEExplicitSolidSolver2::CheckNodalState() const
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    int badNode = mesh.Nodes();

#pragma omp parallel for reduction(min:badNode)
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        const vec3d u  = node.get_vec3d(m_dofU[0],  m_dofU[1],  m_dofU[2]);
        const vec3d v  = node.get_vec3d(m_dofV[0],  m_dofV[1],  m_dofV[2]);
        const vec3d su = node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]);
        const vec3d sv = node.get_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2]);
        const vec3d q  = node.get_vec3d(m_dofQ[0],  m_dofQ[1],  m_dofQ[2]);
        const vec3d w  = node.get_vec3d(m_dofBW[0], m_dofBW[1], m_dofBW[2]);
        const vec3d sa = node.get_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2]);
        const vec3d ba = node.get_vec3d(m_dofBA[0], m_dofBA[1], m_dofBA[2]);
        if (!IsFiniteVec3d(u) || !IsFiniteVec3d(v) || !IsFiniteVec3d(node.m_at) ||
            !IsFiniteVec3d(su) || !IsFiniteVec3d(sv) || !IsFiniteVec3d(sa) ||
            !IsFiniteVec3d(q) || !IsFiniteVec3d(w) || !IsFiniteVec3d(ba))
            badNode = std::min(badNode, i);
    }

    if (badNode < mesh.Nodes())
    {
        feLogError("Explicit solid update produced a non-finite state at node %d.", badNode + 1);
        return false;
    }
    return true;
}

void FEExplicitSolidSolver2::LogState(const char* label) const
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    double umax = 0.0, vmax = 0.0, amax = 0.0;

#pragma omp parallel for reduction(max:umax,vmax,amax)
    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        FENode& node = mesh.Node(i);
        umax = std::max(umax, node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]).norm());
        umax = std::max(umax, node.get_vec3d(m_dofSU[0], m_dofSU[1], m_dofSU[2]).norm());
        umax = std::max(umax, node.get_vec3d(m_dofQ[0], m_dofQ[1], m_dofQ[2]).norm());
        vmax = std::max(vmax, node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]).norm());
        vmax = std::max(vmax, node.get_vec3d(m_dofSV[0], m_dofSV[1], m_dofSV[2]).norm());
        vmax = std::max(vmax, node.get_vec3d(m_dofBW[0], m_dofBW[1], m_dofBW[2]).norm());
        amax = std::max(amax, node.m_at.norm());
        amax = std::max(amax, node.get_vec3d(m_dofSA[0], m_dofSA[1], m_dofSA[2]).norm());
        amax = std::max(amax, node.get_vec3d(m_dofBA[0], m_dofBA[1], m_dofBA[2]).norm());
    }

    feLog("%s: max primary=%.16g, max rate=%.16g, max acceleration=%.16g\n",
        label, umax, vmax, amax);
}

bool FEExplicitSolidSolver2::Quasin()
{
    FEModel& fem = *GetFEModel();
    const double dt = fem.GetTime().timeIncrement;

    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        feLogError("The explicit solid solver requires a positive finite time increment.");
        return false;
    }

    m_niter = 0; m_nrhs = 0; m_nref = 0; m_ntotref = 0;

    PrepStep();
    const vector<double> prescribed = m_ui;

    CaptureAcceptedVelocity();
    InsertPrescribedPrimaryDofs(prescribed, m_velocity, 1.0 / dt);

    vector<double> zeroU(m_neq, 0.0), a0;
    if (!EvaluateAcceleration(zeroU, m_velocity, a0)) return false;

    vector<double> du(m_neq, 0.0), vPredict(m_neq, 0.0);
    double maxDu = 0.0;
#pragma omp parallel for reduction(max:maxDu)
    for (int i = 0; i < m_neq; ++i)
    {
        if (m_lumpedMass[i] > m_massFloor)
        {
            du[i] = dt * m_velocity[i] + 0.5 * dt * dt * a0[i];
            vPredict[i] = m_velocity[i] + dt * a0[i];
            maxDu = std::max(maxDu, std::fabs(du[i]));
        }
    }

    InsertPrescribedPrimaryDofs(prescribed, du, 1.0);
    InsertPrescribedPrimaryDofs(prescribed, vPredict, 1.0 / dt);

    vector<double> a1;
    if (!EvaluateAcceleration(du, vPredict, a1)) return false;

    vector<double> vNew(m_neq, 0.0);
    double maxDv = 0.0;
#pragma omp parallel for reduction(max:maxDv)
    for (int i = 0; i < m_neq; ++i)
    {
        if (m_lumpedMass[i] > m_massFloor)
        {
            const double dv = 0.5 * dt * (a0[i] + a1[i]);
            vNew[i] = m_velocity[i] + dv;
            maxDv = std::max(maxDv, std::fabs(dv));
        }
    }
    InsertPrescribedPrimaryDofs(prescribed, vNew, 1.0 / dt);

    m_stageVelocity = vNew;
    m_stageAcceleration = a1;
    vector<double> acceptedIncrement(du);
    Update(acceptedIncrement);

    double minJ = 0.0;
    if (!CheckElementJacobians("accepted", minJ)) return false;
    if (!CheckNodalState()) return false;

    UpdateIncrements(m_Ui, acceptedIncrement, false);
    UpdateIncrements(m_Ut, m_Ui, true);
    zero(m_Ui);

    m_velocity = vNew;
    m_acceleration = a1;

    feLog("Explicit solid/shell/beam update: max|dq|=%.16g, max|dv|=%.16g, minJ=%.16g\n",
        maxDu, maxDv, minJ);
    LogState("Explicit solid post-update state");

    m_niter = 1;
    fem.DoCallback(CB_MINOR_ITERS);
    fem.DoCallback(CB_QUASIN_CONVERGED);
    return true;
}

void FEExplicitSolidSolver2::Serialize(DumpStream& ar)
{
    FESolidSolver2::Serialize(ar);
    ar & m_lumpedMass & m_velocity & m_acceleration;
    ar & m_stageVelocity & m_stageAcceleration;
    ar & m_dynDamping & m_massFloor & m_minJ;
}
