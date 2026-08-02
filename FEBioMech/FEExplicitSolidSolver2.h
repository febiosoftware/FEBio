#pragma once

#include "FESolidSolver2.h"
#include "febiomech_api.h"

// Explicit velocity-Verlet / central-difference solver for three-dimensional
// deformable solid domains.
//
// The solver advances
//
//     u_dot = v,
//     M_L v_dot = R(u, v),
//
// where M_L is the row-sum-lumped translational mass. It derives from
// FESolidSolver2 so that FEBio's equation numbering, prescribed-displacement
// handling, load preparation, domain updates, and accepted-step bookkeeping
// remain available.
//
// Current scope:
//   * dynamic analyses only;
//   * FEElasticSolidDomain domains only;
//   * translational displacement equations only;
//   * no shells, beams, rigid bodies, nonlinear constraints, or contact;
//   * no mass scaling or automatic critical-time-step estimator yet.
//
// OpenMP follows FEBio's pragma-based convention.
class FEBIOMECH_API FEExplicitSolidSolver2 : public FESolidSolver2
{
public:
    explicit FEExplicitSolidSolver2(FEModel* fem);
    ~FEExplicitSolidSolver2() override = default;

    bool Init() override;
    void PrepStep() override;
    bool Quasin() override;
    void UpdateKinematics(vector<double>& ui) override;
    void Serialize(DumpStream& ar) override;

protected:
    bool BuildLumpedMass();

    bool EvaluateAcceleration(const vector<double>& stageU,
        const vector<double>& stageV,
        vector<double>& acceleration);

    bool FormAcceleration(const vector<double>& R,
        vector<double>& acceleration) const;

    void CaptureAcceptedVelocity();
    void InsertPrescribedDisplacement(const vector<double>& prescribed,
        vector<double>& increment) const;
    void InsertPrescribedVelocity(const vector<double>& prescribed,
        vector<double>& velocity, double dt) const;

    bool CheckNodalState() const;
    bool CheckElementJacobians(const char* stage, double& minJ) const;
    void LogState(const char* label) const;

protected:
    vector<double> m_lumpedMass;
    vector<double> m_velocity;
    vector<double> m_acceleration;

    // State used by the overridden UpdateKinematics during a trial stage.
    vector<double> m_stageVelocity;
    vector<double> m_stageAcceleration;

public:
    double m_dynDamping;
    double m_massFloor;
    double m_minJ;

    DECLARE_FECORE_CLASS();
};
