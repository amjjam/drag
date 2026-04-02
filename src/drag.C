#include "drag.H"
#include <cmath>

// =============================
// Time Management (Epoch)
// =============================

Epoch::Epoch() : jd_utc_(0.0) {}

Epoch::Epoch(double jd_utc) : jd_utc_(jd_utc) {}

// Inital convertion to UTC Julian Date
Epoch::Epoch(const std::tm& tm) {
    int year = tm.tm_year + 1900;
    int month = tm.tm_mon + 1;
    int day = tm.tm_mday;
    int hour = tm.tm_hour;
    int min = tm.tm_min;
    int sec = tm.tm_sec;

    if (month <= 2) {
        year -= 1;
        month += 12;
    }

    double A = std::floor(year / 100.0);
    double B = 2 - A + std::floor(A / 4.0);

    double jd = std::floor(365.25 * (year + 4716)) +
                std::floor(30.6001 * (month + 1)) +
                day + B - 1524.5;

    double fraction = (hour + min / 60.0 + sec / 3600.0) / 24.0;
    jd_utc_ = jd + fraction;
}

// Operator Overload
Epoch Epoch::operator+(double seconds) const {
    // 86400 seconds in a day. Add the fraction of a day to JD.
    return Epoch(jd_utc_ + (seconds / 86400.0));
}

double Epoch::toJD() const {
    return jd_utc_;
}

// Convert UTC Julian Date to Terrestrial Time Modified Julian Date
double Epoch::toMjdTT() const {
    double mjd_utc = jd_utc_ - 2400000.5;
    // TT is roughly UTC + 69.184 seconds (37s leap seconds + 32.184s offset)
    return mjd_utc + (69.184 / 86400.0); 
}

// Get GMST anlge
double Epoch::toGMST() const {
    double T_UT1 = (jd_utc_ - 2451545.0) / 36525.0;
    double gmst_seconds = 67310.54841 + 
                          (876600.0 * 3600 + 8640184.812866) * T_UT1 + 
                          0.093104 * T_UT1 * T_UT1 - 
                          6.2e-6 * T_UT1 * T_UT1 * T_UT1;

    gmst_seconds = std::fmod(gmst_seconds, 86400.0);
    if (gmst_seconds < 0) gmst_seconds += 86400.0;

    return gmst_seconds * (2.0 * M_PI / 86400.0);
}

std::time_t Epoch::toTimeT() const {
    // JD to Unix time (seconds since 1970-01-01)
    // JD 2440587.5 is 1970-01-01 00:00:00 UTC
    return static_cast<std::time_t>(std::round((jd_utc_ - 2440587.5) * 86400.0));
}

// =============================
// GravityAccel
// =============================

GravityAccel::GravityAccel(double mu) : mu_(mu) {}

Eigen::Vector3d GravityAccel::computeAcceleration(
    const Spacecraft&,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    Epoch
) const {
    double r = pos.norm();
    return -mu_ * pos / (r * r * r);
}

// ==============================
// EGM2008GravityAccel
// ==============================

EGM2008GravityAccel::EGM2008GravityAccel(const std::string& model_name, int max_degree, int max_order)
    	: grav_(model_name, "", max_degree, max_order) {}
Eigen::Vector3d EGM2008GravityAccel::computeAcceleration(
	const Spacecraft&,
	const Eigen::Vector3d& pos_eci,
	const Eigen::Vector3d&,
	Epoch t
) const {
	const double rotation_angle = t.toGMST();
	Eigen::Matrix3d C_e_i;
	C_e_i <<  cos(rotation_angle), sin(rotation_angle), 0,
	         -sin(rotation_angle), cos(rotation_angle), 0,
	          0,                   0,                 1;

	Eigen::Vector3d pos_ecef = C_e_i * pos_eci;
	Eigen::Vector3d g_ecef;
	grav_.V(pos_ecef.x(), pos_ecef.y(), pos_ecef.z(), g_ecef.x(), g_ecef.y(), g_ecef.z());
	Eigen::Matrix3d C_i_e = C_e_i.transpose();
	Eigen::Vector3d g_eci = C_i_e * g_ecef;
	return g_eci;
}

// =============================
// SunGravityAccel
// =============================

SunGravityAccel::SunGravityAccel() {}

Eigen::Vector3d SunGravityAccel::computeAcceleration(
	const Spacecraft&,
	const Eigen::Vector3d& pos_eci,
	const Eigen::Vector3d&,
	Epoch t
) const {
	std::time_t current_time = t.toTimeT();
	Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);
	Eigen::Vector3d pos_wrt_sun = pos_eci - sun_pos;

	return (-1.32712438e+20) * (pos_wrt_sun / pow(pos_wrt_sun.norm(),3) + sun_pos / pow(sun_pos.norm(),3));
}

// =============================
// MoonGravityAccel
// =============================
/*
MoonGravityAccel::MoonGravityAccel(std::time_t epoch) : epoch_(epoch) {}

Eigen::Vector3d MoonGravityAccel::computeAcceleration(
        const Spacecraft&,
        const Eigen::Vector3d& pos_eci,
        const Eigen::Vecotr3d&,
        const t
) const {
        std::time_t current_time = epoch_ + static_cast<long>(t);
        Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);
        Eigen::Vector3d pos_wrt_sun = pos_eci - sun_pos;

        return (-GM) * (pos_wrt_sun / pow(pos_wrt_sun.norm(),3) + sun_pos / pow(sun_pos.norm(),3));
}
*/

// =============================
// DragAccel
// =============================

DragAccel::DragAccel(double rho0) : rho0_(rho0) {}

Eigen::Vector3d DragAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d&,
    const Eigen::Vector3d& vel,
    Epoch
) const {
    double v = vel.norm();
    if (v == 0.0) return Eigen::Vector3d::Zero();
    return -0.5 * sc.Cd() * sc.area() * rho0_ * v / sc.mass() * vel;
}


// =============================
// MSISDragAccel
// =============================

MsisDragAccel::MsisDragAccel(double f107a, double f107, double ap)
    : f107_avg(f107a),
      f107_daily(f107),
      ap_val(ap),
      earth(GeographicLib::Constants::WGS84_a(), GeographicLib::Constants::WGS84_f()),
      omega(GeographicLib::Constants::WGS84_omega())
{
    // --- Initialize the MSIS model ---
    // This constructor requires a flags array.
    std::array<int, 24> flags;
    
    // Flag 0 = 1 means use SI units (kg/m^3), which is what we want
    flags[0] = 1; 
    
    // Set all other flags to 'on' (standard practice)
    for (int i = 1; i < 24; i++) {
        flags[i] = 1;
    }

    // Create the model object and store it in our smart pointer
    msis_model_ptr = std::make_unique< ::atmos::CNrlmsise00 >(flags);

    // --- Initialize the AP array ---
    // Fill the entire array with the single daily 'ap' value
    ap_array.fill(ap_val);
}

Eigen::Vector3d MsisDragAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos_eci,
    const Eigen::Vector3d& vel_eci,
    Epoch t
) const {
    
    // === 1. Calculate Current Time ===
    std::time_t current_time_sec = t.toTimeT();

    std::tm current_time_utc;
    gmtime_r(&current_time_sec, &current_time_utc); // Thread-safe UTC time

    int doy = current_time_utc.tm_yday + 1; // Day of year (tm_yday is 0-365)
    double sec_of_day = current_time_utc.tm_hour * 3600.0 + 
                       current_time_utc.tm_min * 60.0 + 
                       current_time_utc.tm_sec;

    // === 2. Convert ECI Position to Geodetic (Lat, Lon, Alt) ===
    double rotation_angle = t.toGMST(); // Simplified rotation
    Eigen::Matrix3d C_e_i;  // ECI to ECEF
    C_e_i <<  cos(rotation_angle), sin(rotation_angle), 0,
             -sin(rotation_angle), cos(rotation_angle), 0,
              0,                   0,                   1;
    Eigen::Vector3d r_ecef = C_e_i * pos_eci;

    double lat, lon, alt_m;
    earth.Reverse(r_ecef.x(), r_ecef.y(), r_ecef.z(), lat, lon, alt_m);
    double alt_km = alt_m / 1000.0; // Library requires altitude in [km]

    // === 3. Call MSIS ===
    // Use the simple 'density' function from the library
    double rho = msis_model_ptr->density(
        doy,
        sec_of_day,
        alt_km,
        lat,
        lon,
        f107_avg,
        f107_daily,
        ap_array // Pass the stored AP array
    );
    // 'rho' is now in [kg/m^3] because we set flag 0 = 1.

    // === 4. Calculate Drag ===
    Eigen::Vector3d v_atm_eci = Eigen::Vector3d(-omega * pos_eci.y(), omega * pos_eci.x(), 0.0);
    Eigen::Vector3d v_rel = vel_eci - v_atm_eci;
    double v_mag = v_rel.norm();

    if (v_mag == 0.0) return Eigen::Vector3d::Zero();

    double B = -0.5 * (sc.Cd() * sc.area() / sc.mass()) * rho * v_mag;
    Eigen::Vector3d a_drag_eci = B * v_rel;

    return a_drag_eci;
}


// =============================
// SolarRadiationAccel
// =============================

SolarRadiationAccel::SolarRadiationAccel(double P0) : P0_(P0) {}

Eigen::Vector3d SolarRadiationAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d&,
    const Eigen::Vector3d&,
    Epoch
) const {
    Eigen::Vector3d sunDir(1.0, 0.0, 0.0); // assume Sun in +x
    return (P0_ * sc.Cr() * sc.area() / sc.mass()) * sunDir;
}

// =============================
// DynamicSolarAccel
// =============================

DynamicSolarAccel::DynamicSolarAccel(double P0, double earth_radius)
    : P0_(P0), earth_radius_(earth_radius) {}

Eigen::Vector3d DynamicSolarAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    Epoch t
) const {
    std::time_t current_time = t.toTimeT();
    Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);
    Eigen::Vector3d sun_dir = sun_pos.normalized(); // Unit vector pointing to Sun

    double s = pos.dot(sun_dir);
    if (s < 0) {
	Eigen::Vector3d dist_vec = pos - (s * sun_dir);
	double dist_sq = dist_vec.squaredNorm();
 	if (dist_sq < (earth_radius_ * earth_radius_)) {
	    return Eigen::Vector3d::Zero();
        }
    }
    return (P0_ * sc.Cr() * sc.area() / sc.mass()) * (-sun_dir);
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
    Epoch t
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

State Dynamics::operator()(const State& s, Epoch t) const {
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
                          Epoch t,
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
                                              Epoch t0,
                                              double dt,
                                              int steps) const {
    std::vector<State> trajectory;
    trajectory.reserve(steps+1);
    trajectory.push_back(s0);

    State s = s0;
    Epoch t = t0;
    for (int i = 0; i < steps; ++i) {
        s = integrator_->step(dyn_, s, t, dt);
        t = t + dt;
        trajectory.push_back(s);
    }
    return trajectory;
}
