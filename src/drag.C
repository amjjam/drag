#include "drag.H"
#include <cmath>

// =============================
// GravityAccel
// =============================

GravityAccel::GravityAccel(double mu) : mu_(mu) {}

Eigen::Vector3d GravityAccel::computeAcceleration(
    const Spacecraft&,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    double
) const {
    double r = pos.norm();
    return -mu_ * pos / (r * r * r);
}

// =============================
// DragAccel
// =============================

DragAccel::DragAccel(double rho0) : rho0_(rho0) {}

Eigen::Vector3d DragAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d&,
    const Eigen::Vector3d& vel,
    double
) const {
    double v = vel.norm();
    if (v == 0.0) return Eigen::Vector3d::Zero();
    return -0.5 * sc.Cd() * sc.area() * rho0_ * v / sc.mass() * vel;
}

// =============================
// SolarRadiationAccel
// =============================

SolarRadiationAccel::SolarRadiationAccel(double P0) : P0_(P0) {}

Eigen::Vector3d SolarRadiationAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d&,
    const Eigen::Vector3d&,
    double
) const {
    Eigen::Vector3d sunDir(1.0, 0.0, 0.0); // assume Sun in +x
    return (P0_ * sc.Cr() * sc.area() / sc.mass()) * sunDir;
}

// =============================
// AccelAggregator
// =============================

void AccelAggregator::addModel(std::shared_ptr<AccelModel> model) {
    models_.push_back(model);
}

Eigen::Vector3d AccelAggregator::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d& vel,
    double t
) const {
    Eigen::Vector3d total = Eigen::Vector3d::Zero();
    for (const auto& m : models_) {
        total += m->computeAcceleration(sc, pos, vel, t);
    }
    return total;
}

// =============================
// Dynamics
// =============================

Dynamics::Dynamics(const Spacecraft& sc, const AccelAggregator& accels)
    : sc_(sc), accels_(accels) {}

State Dynamics::operator()(const State& s, double t) const {
    State deriv;
    deriv.pos = s.vel;
    deriv.vel = accels_.computeAcceleration(sc_, s.pos, s.vel, t);
    return deriv;
}

// =============================
// RK4Integrator
// =============================

State RK4Integrator::step(const Dynamics& dyn,
                          const State& s,
                          double t,
                          double dt) const {
    State k1 = dyn(s, t);
    State k2 = dyn({s.pos + 0.5 * dt * k1.pos,
                    s.vel + 0.5 * dt * k1.vel}, t + 0.5 * dt);
    State k3 = dyn({s.pos + 0.5 * dt * k2.pos,
                    s.vel + 0.5 * dt * k2.vel}, t + 0.5 * dt);
    State k4 = dyn({s.pos + dt * k3.pos,
                    s.vel + dt * k3.vel}, t + dt);

    State out;
    out.pos = s.pos + (dt/6.0) * (k1.pos + 2.0*k2.pos + 2.0*k3.pos + k4.pos);
    out.vel = s.vel + (dt/6.0) * (k1.vel + 2.0*k2.vel + 2.0*k3.vel + k4.vel);
    return out;
}

// =============================
// OrbitPropagator
// =============================

OrbitPropagator::OrbitPropagator(const Spacecraft& sc,
                                 const AccelAggregator& accels,
                                 std::shared_ptr<Integrator> integrator)
    : dyn_(sc, accels), integrator_(std::move(integrator)) {}

std::vector<State> OrbitPropagator::propagate(const State& s0,
                                              double t0,
                                              double dt,
                                              int steps) const {
    std::vector<State> trajectory;
    trajectory.reserve(steps+1);
    trajectory.push_back(s0);

    State s = s0;
    double t = t0;
    for (int i = 0; i < steps; ++i) {
        s = integrator_->step(dyn_, s, t, dt);
        t += dt;
        trajectory.push_back(s);
    }
    return trajectory;
}
