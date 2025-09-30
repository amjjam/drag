#ifndef _DRAG_H_
#define _DRAG_H_

#include <Eigen/Dense>
#include <vector>
#include <memory>

// =============================
// Data Structures
// =============================

struct State {
    Eigen::Vector3d pos; // position [m]
    Eigen::Vector3d vel; // velocity [m/s]
};

struct Spacecraft {
    double mass;  // [kg]
    double area;  // [m^2]
    double Cd;    // drag coefficient
    double Cr;    // reflectivity coefficient
};

// =============================
// Acceleration Models
// =============================

class AccelModel {
public:
    virtual ~AccelModel() = default;
    virtual Eigen::Vector3d computeAcceleration(
        const Spacecraft& sc,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double t
    ) const = 0;
};

class GravityAccel : public AccelModel {
public:
    explicit GravityAccel(double mu);
    Eigen::Vector3d computeAcceleration(
        const Spacecraft& sc,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double t
    ) const override;
private:
    double mu_;
};

class DragAccel : public AccelModel {
public:
    explicit DragAccel(double rho0);
    Eigen::Vector3d computeAcceleration(
        const Spacecraft& sc,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double t
    ) const override;
private:
    double rho0_;
};

class SolarRadiationAccel : public AccelModel {
public:
    explicit SolarRadiationAccel(double P0);
    Eigen::Vector3d computeAcceleration(
        const Spacecraft& sc,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double t
    ) const override;
private:
    double P0_;
};

class AccelAggregator : public AccelModel {
public:
    void addModel(std::shared_ptr<AccelModel> model);
    Eigen::Vector3d computeAcceleration(
        const Spacecraft& sc,
        const Eigen::Vector3d& pos,
        const Eigen::Vector3d& vel,
        double t
    ) const override;
private:
    std::vector<std::shared_ptr<AccelModel>> models_;
};

// =============================
// Dynamics
// =============================

class Dynamics {
public:
    Dynamics(const Spacecraft& sc, const AccelAggregator& accels);
    State operator()(const State& s, double t) const;
private:
    const Spacecraft& sc_;
    const AccelAggregator& accels_;
};

// =============================
// Integrators
// =============================

class Integrator {
public:
    virtual ~Integrator() = default;
    virtual State step(const Dynamics& dyn,
                       const State& s,
                       double t,
                       double dt) const = 0;
};

class RK4Integrator : public Integrator {
public:
    State step(const Dynamics& dyn,
               const State& s,
               double t,
               double dt) const override;
};

// =============================
// Propagator
// =============================

class OrbitPropagator {
public:
    OrbitPropagator(const Spacecraft& sc,
                    const AccelAggregator& accels,
                    std::shared_ptr<Integrator> integrator);

    std::vector<State> propagate(const State& s0,
                                 double t0,
                                 double dt,
                                 int steps) const;
private:
    Dynamics dyn_;
    std::shared_ptr<Integrator> integrator_;
};

#endif /* _DRAG_H_ */

