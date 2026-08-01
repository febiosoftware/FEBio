#pragma once

#include "FEFluidSolver.h"
#include "febiofluid_api.h"

// Explicit classical fourth-order Runge-Kutta solver for FEBio's
// velocity-dilatation fluid formulation.
//
// The class deliberately derives from FEFluidSolver so that equation
// registration, prescribed-DOF handling, global solution vectors,
// kinematic updates, constraints, loads, and residual assembly are inherited
// from the standard fluid solver. Only the nonlinear Newton iteration is
// replaced by a four-stage explicit RK4 update using row-sum-lumped diagonal
// operators.
//
// OpenMP is used through the same pragma-based approach used elsewhere in
// FEBio. No _OPENMP checks, runtime thread-count calls, or solver-specific
// OpenMP parameters are required. Compilers without OpenMP ignore the pragmas.
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

    // Evaluate y_dot = M_L^{-1} R at the state m_Ut + stageIncrement.
    bool EvaluateStage(const vector<double>& stageIncrement, vector<double>& rate);

    // Convert the assembled residual to a rate vector. Prescribed and other
    // equations without a positive primary-fluid capacity receive zero rate.
    bool FormExplicitRate(const vector<double>& R, vector<double>& rate) const;

    // Insert the standard FEBio-prescribed increments into an increment vector.
    // scale=0.5 is used at RK4 midpoint stages and scale=1 at the final stage.
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

public:
    double m_dynDamping;    // multiplies all explicit rates; normally 1
    double m_capacityFloor; // smallest accepted positive lumped coefficient

    DECLARE_FECORE_CLASS();
};
