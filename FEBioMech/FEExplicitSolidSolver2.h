#pragma once

#include "FESolidSolver2.h"
#include "febiomech_api.h"

// Explicit velocity-Verlet / central-difference solver for deformable solids,
// shells, and beams. Penalty contact is supported through the ordinary
// FESolidSolver2 residual and contact lifecycle. Rigid bodies, nonlinear
// constraint equations, and augmented-Lagrangian contact are not supported.
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
    void InsertPrescribedPrimaryDofs(const vector<double>& prescribed,
        vector<double>& values, double velocityScale) const;

    void GatherNodeRate(const FENode& node, const FEDofList& primary,
        const FEDofList& rate, vector<double>& globalRate) const;
    void SetStageNodeRate(FENode& node, const FEDofList& primary,
        const FEDofList& rate, const FEDofList* accelerationRate) const;

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
