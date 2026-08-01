#pragma once

#include "FEFluidSolver.h"
#include "febiofluid_api.h"

// Explicit classical fourth-order Runge-Kutta solver for FEBio's
// velocity-dilatation fluid formulation.
//
// For robustness under finite or large changes in fluid volume ratio, the
// solver integrates q = ln(J), J = 1 + e, rather than integrating e additively.
// The FEBio primary variable remains fluid dilatation e, so each RK stage is
// mapped back exactly with
//
//     J_stage = J_n exp(Delta q_stage),
//     e_stage = J_stage - 1.
//
// This multiplicative update preserves J > 0 whenever the logarithmic stage
// increment is finite. Velocity remains additively integrated.
//
// The class derives from FEFluidSolver so equation registration, prescribed
// DOFs, global solution vectors, constraints, loads, and residual assembly are
// inherited from the standard fluid solver.
class FEBIOFLUID_API FEExplicitFluidSolver : public FEFluidSolver
{
public:
    explicit FEExplicitFluidSolver(FEModel* fem);
    ~FEExplicitFluidSolver() override = default;

    bool Init() override;
    bool Quasin() override;
    void Serialize(DumpStream& ar) override;

protected:
    bool BuildLumpedOperators();

    // Evaluate the mixed rate vector at m_Ut + stageIncrement.
    // Velocity entries contain v_dot. Dilatation entries are converted from
    // e_dot to q_dot = e_dot/J at the current RK stage.
    bool EvaluateStage(const vector<double>& stageIncrement, vector<double>& rate);

    // Convert the assembled residual to rates. Before the logarithmic
    // transformation, dilatation entries contain e_dot.
    bool FormExplicitRate(const vector<double>& R, vector<double>& rate) const;

    // Convert free dilatation equation rates from e_dot to q_dot=e_dot/J
    // using the current nodal stage state.
    bool ConvertDilatationRatesToLogRates(vector<double>& rate) const;

    // Capture J_n at the start of the accepted time step.
    bool CaptureBaseVolumeRatios();

    // Form an RK stage increment from a mixed rate vector. Velocity is updated
    // additively. Free dilatation is updated multiplicatively through ln(J).
    bool BuildStageIncrement(const vector<double>& rate, double factor,
        const vector<double>& prescribed, double prescribedScale,
        vector<double>& increment) const;

    // Form the accepted RK4 increment using the weighted velocity rates and
    // weighted logarithmic volume-ratio rates.
    bool BuildFinalIncrement(const vector<double>& k1, const vector<double>& k2,
        const vector<double>& k3, const vector<double>& k4, double dt,
        const vector<double>& prescribed, vector<double>& increment,
        double& maxRate, double& maxIncrement) const;

    // Insert standard FEBio-prescribed increments. These remain controlled by
    // the prescribed boundary-condition implementation.
    void InsertPrescribedIncrements(const vector<double>& prescribed,
        vector<double>& increment, double scale) const;

    void ZeroPrimaryRates();
    bool CheckExplicitState() const;
    void LogPrimaryState(const char* label) const;

protected:
    // Row-sum-lumped diagonal coefficient indexed by global equation number.
    // Velocity entries are integral(rho*N_a dv); dilatation entries are
    // integral(N_a/J dv). Only free primary-fluid equations are populated.
    vector<double> m_capacity;

    // Nodal volume ratio at the beginning of the current accepted step.
    vector<double> m_baseJ;

public:
    double m_dynDamping;       // multiplies all explicit rates; normally 1
    double m_capacityFloor;    // smallest accepted positive lumped coefficient
    double m_minVolumeRatio;   // minimum accepted J at stages and final state
    double m_maxLogIncrement;  // maximum |Delta ln J| per RK stage; <=0 disables

    DECLARE_FECORE_CLASS();
};
