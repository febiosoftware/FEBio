#include "stdafx.h"
#include "FEExplicitSolidSolver2.h"
#include "FEElasticSolidDomain.h"
#include "FEElasticMaterial.h"
#include "FESolidMaterial.h"
#include "FESolidAnalysis.h"

#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/FENode.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FEBoundaryCondition.h>
#include <FECore/FELinearConstraintManager.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FEModelLoad.h>
#include <FECore/log.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
    inline bool IsFiniteVec3d(const vec3d& v)
    {
        return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
    }
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

    // Initial accelerations are calculated from the explicit residual and the
    // lumped mass, not by the consistent-mass linear solve in FESolidSolver2.
    m_init_accelerations = false;

    // This solver never reforms a stiffness matrix during time stepping.
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
        feLogError("The explicit solid solver does not currently support nonlinear-constraint equations.");
        return false;
    }

    if (fem.SurfacePairConstraints() != 0)
    {
        feLogError("The explicit solid solver does not currently support contact or surface-pair constraints.");
        return false;
    }

    for (int i = 0; i < mesh.Nodes(); ++i)
    {
        if (mesh.Node(i).m_rid >= 0)
        {
            feLogError("The explicit solid solver does not currently support rigid bodies or rigid nodes.");
            return false;
        }
    }

    for (int i = 0; i < mesh.Domains(); ++i)
    {
        FEElasticSolidDomain* dom = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(i));
        if (dom == nullptr)
        {
            feLogError("The explicit solid solver currently supports FEElasticSolidDomain only (domain %d is unsupported).", i + 1);
            return false;
        }
        dom->SetDynamicUpdateFlag(true);
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

    // FEBio stores prescribed displacement increments in m_ui.
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

    // The accepted state is already present on the nodes. Update material-point
    // history before evaluating the first explicit residual of this step.
    UpdateModel();

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

    m_lumpedMass.assign(m_neq, 0.0);

    for (int nd = 0; nd < mesh.Domains(); ++nd)
    {
        FEElasticSolidDomain* dom = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
        if (dom == nullptr) return false;

        FESolidMaterial* mat = dynamic_cast<FESolidMaterial*>(dom->GetMaterial());
        if (mat == nullptr)
        {
            feLogError("Solid domain %d does not reference an FESolidMaterial.", nd + 1);
            return false;
        }

        const int NE = dom->Elements();
        int failed = 0;

#pragma omp parallel for schedule(static) reduction(|:failed)
        for (int iel = 0; iel < NE; ++iel)
        {
            FESolidElement& el = dom->Element(iel);
            if (!el.isActive()) continue;

            const int neln = el.Nodes();
            vector<double> me(neln, 0.0);
            bool elementOK = true;

            for (int n = 0; n < el.GaussPoints(); ++n)
            {
                FEMaterialPoint& mp = *el.GetMaterialPoint(n);
                const double rho = mat->Density(mp);
                const double dV0 = dom->detJ0(el, n) * el.GaussWeights()[n];
                const double* H = el.H(n);

                if (!(rho > 0.0) || !std::isfinite(rho) ||
                    !(dV0 > 0.0) || !std::isfinite(dV0))
                {
                    elementOK = false;
                    break;
                }

                // Row sum of integral rho N_a N_b dV0.
                for (int a = 0; a < neln; ++a)
                    me[a] += rho * H[a] * dV0;
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
                    const int eq = node.m_ID[m_dofU[k]];
                    if ((eq >= 0) && (eq < m_neq))
                    {
#pragma omp atomic
                        m_lumpedMass[eq] += me[a];
                    }
                }
            }
        }

        if (failed)
        {
            feLogError("Invalid density or reference Jacobian while assembling lumped mass in domain %d.", nd + 1);
            return false;
        }
    }

    int active = 0;
    double mmin = std::numeric_limits<double>::max();
    double mmax = 0.0;
#pragma omp parallel for reduction(+:active) reduction(min:mmin) reduction(max:mmax)
    for (int i = 0; i < m_neq; ++i)
    {
        const double m = m_lumpedMass[i];
        if (m > m_massFloor)
        {
            active += 1;
            mmin = std::min(mmin, m);
            mmax = std::max(mmax, m);
        }
    }

    if (active == 0)
    {
        feLogError("No positive translational lumped masses were assembled.");
        return false;
    }

    feLog("Explicit solid lumped mass: active=%d, min=%.16g, max=%.16g\n", active, mmin, mmax);
    return true;
}

void FEExplicitSolidSolver2::CaptureAcceptedVelocity()
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    std::fill(m_velocity.begin(), m_velocity.end(), 0.0);

    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        const vec3d v = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
        for (int k = 0; k < 3; ++k)
        {
            const int id = node.m_ID[m_dofU[k]];
            const int eq = (id < -1 ? -id - 2 : id);
            if ((eq >= 0) && (eq < m_neq))
                m_velocity[eq] = (k == 0 ? v.x : (k == 1 ? v.y : v.z));
        }
    }
}

void FEExplicitSolidSolver2::InsertPrescribedDisplacement(
    const vector<double>& prescribed, vector<double>& increment) const
{
    FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        for (int k = 0; k < 3; ++k)
        {
            const int id = node.m_ID[m_dofU[k]];
            if (id < -1)
            {
                const int eq = -id - 2;
                if ((eq >= 0) && (eq < (int)increment.size()))
                    increment[eq] = prescribed[eq];
            }
        }
    }
}

void FEExplicitSolidSolver2::InsertPrescribedVelocity(
    const vector<double>& prescribed, vector<double>& velocity, double dt) const
{
    if (!(dt > 0.0)) return;
    FEMesh& mesh = GetFEModel()->GetMesh();
    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        for (int k = 0; k < 3; ++k)
        {
            const int id = node.m_ID[m_dofU[k]];
            if (id < -1)
            {
                const int eq = -id - 2;
                if ((eq >= 0) && (eq < (int)velocity.size()))
                    velocity[eq] = prescribed[eq] / dt;
            }
        }
    }
}

void FEExplicitSolidSolver2::UpdateKinematics(vector<double>& ui)
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();

    vector<double> U(m_Ut.size(), 0.0);
    const int neq = (int)U.size();
#pragma omp parallel for
    for (int i = 0; i < neq; ++i)
        U[i] = ui[i] + m_Ui[i] + m_Ut[i];

    scatter3(U, mesh, m_dofU[0], m_dofU[1], m_dofU[2]);

    for (int i = 0; i < fem.BoundaryConditions(); ++i)
    {
        FEBoundaryCondition& bc = *fem.BoundaryCondition(i);
        if (bc.IsActive()) bc.Update();
    }

    FELinearConstraintManager& LCM = fem.GetLinearConstraintManager();
    if (LCM.LinearConstraints() > 0) LCM.Update();

    const int NN = mesh.Nodes();
#pragma omp parallel for
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        node.m_rt = node.m_r0 + node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);

        vec3d v(0, 0, 0), a(0, 0, 0);
        for (int k = 0; k < 3; ++k)
        {
            const int id = node.m_ID[m_dofU[k]];
            const int eq = (id < -1 ? -id - 2 : id);
            if ((eq >= 0) && (eq < (int)m_stageVelocity.size()))
            {
                const double vk = m_stageVelocity[eq];
                const double ak = m_stageAcceleration[eq];
                if (k == 0) { v.x = vk; a.x = ak; }
                else if (k == 1) { v.y = vk; a.y = ak; }
                else { v.z = vk; a.z = ak; }
            }
        }
        node.set_vec3d(m_dofV[0], m_dofV[1], m_dofV[2], v);
        node.m_at = a;
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
            if (!std::isfinite(a))
            {
                invalidEq = std::min(invalidEq, i);
                continue;
            }
            acceleration[i] = a;
            active += 1;
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
        FEElasticSolidDomain* dom = dynamic_cast<FEElasticSolidDomain*>(&mesh.Domain(nd));
        if (dom == nullptr) continue;

        for (int iel = 0; iel < dom->Elements(); ++iel)
        {
            FESolidElement& el = dom->Element(iel);
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
                    badDomain = nd;
                    badElement = iel;
                    badPoint = n;
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
    return true;
}

bool FEExplicitSolidSolver2::EvaluateAcceleration(
    const vector<double>& stageU,
    const vector<double>& stageV,
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
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    const int NN = mesh.Nodes();
    int badNode = NN;

#pragma omp parallel for reduction(min:badNode)
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        const vec3d u = node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]);
        const vec3d v = node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]);
        if (!IsFiniteVec3d(u) || !IsFiniteVec3d(v) || !IsFiniteVec3d(node.m_at))
            badNode = std::min(badNode, i);
    }

    if (badNode < NN)
    {
        feLogError("Explicit solid update produced a non-finite state at node %d.", badNode + 1);
        return false;
    }
    return true;
}

void FEExplicitSolidSolver2::LogState(const char* label) const
{
    FEModel& fem = *GetFEModel();
    FEMesh& mesh = fem.GetMesh();
    const int NN = mesh.Nodes();
    double umax = 0.0, vmax = 0.0, amax = 0.0;

#pragma omp parallel for reduction(max:umax,vmax,amax)
    for (int i = 0; i < NN; ++i)
    {
        FENode& node = mesh.Node(i);
        umax = std::max(umax, node.get_vec3d(m_dofU[0], m_dofU[1], m_dofU[2]).norm());
        vmax = std::max(vmax, node.get_vec3d(m_dofV[0], m_dofV[1], m_dofV[2]).norm());
        amax = std::max(amax, node.m_at.norm());
    }

    feLog("%s: umax=%.16g, vmax=%.16g, amax=%.16g\n", label, umax, vmax, amax);
}

bool FEExplicitSolidSolver2::Quasin()
{
    FEModel& fem = *GetFEModel();
    FETimeInfo& tp = fem.GetTime();
    const double dt = tp.timeIncrement;

    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        feLogError("The explicit solid solver requires a positive finite time increment.");
        return false;
    }

    m_niter = 0;
    m_nrhs = 0;
    m_nref = 0;
    m_ntotref = 0;

    PrepStep();
    const vector<double> prescribed = m_ui;

    // Synchronize the accepted velocity from the state saved by PrepStep.
    CaptureAcceptedVelocity();
    InsertPrescribedVelocity(prescribed, m_velocity, dt);

    vector<double> zeroU(m_neq, 0.0);
    vector<double> a0;
    if (!EvaluateAcceleration(zeroU, m_velocity, a0)) return false;

    vector<double> du(m_neq, 0.0);
    vector<double> vPredict(m_neq, 0.0);
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

    InsertPrescribedDisplacement(prescribed, du);
    InsertPrescribedVelocity(prescribed, vPredict, dt);

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
    InsertPrescribedVelocity(prescribed, vNew, dt);

    // Re-evaluate the accepted state with its accepted velocity and acceleration.
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

    feLog("Explicit solid central-difference update: max|du|=%.16g, max|dv|=%.16g, minJ=%.16g\n",
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
