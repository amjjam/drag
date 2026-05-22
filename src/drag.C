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
	double current_time = t.toJD();
	Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);
	Eigen::Vector3d pos_wrt_sun = pos_eci - sun_pos;

	const double GM_Sun = 1.32712438e+20;

	return (-GM_Sun) * (pos_wrt_sun / pow(pos_wrt_sun.norm(),3) + sun_pos / pow(sun_pos.norm(),3));
}

// =============================
// MoonGravityAccel
// =============================

MoonGravityAccel::MoonGravityAccel() {}

Eigen::Vector3d MoonGravityAccel::computeAcceleration(
        const Spacecraft&,
        const Eigen::Vector3d& pos_eci,
        const Eigen::Vector3d&,
        Epoch t
) const {
        double current_time = t.toMjdTT();
        Eigen::Vector3d moon_pos = MoonEphemeris::getMoonPositionECI(current_time);
        Eigen::Vector3d pos_wrt_moon = pos_eci - moon_pos;

	const double GM_Moon = 398600.4415e+9 / 81.300587;

        return (-GM_Moon) * (pos_wrt_moon / pow(pos_wrt_moon.norm(),3) + moon_pos / pow(moon_pos.norm(),3));
}


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

MsisDragAccel::MsisDragAccel(double f107a, double f107, const std::array<double, 7>& ap_array)
    : f107_avg(f107a),
      f107_daily(f107),
      ap_array_(ap_array),
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
        ap_array_
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
// SRP_CylindricalAccel
// =============================

SRP_CylindricalAccel::SRP_CylindricalAccel(double P0, double earth_radius)
    : P0_(P0), earth_radius_(earth_radius) {}

Eigen::Vector3d SRP_CylindricalAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    Epoch t
) const {
    double current_time = t.toJD();
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
// SRP_CanonicalAccel
// =============================

SRP_CanonicalAccel::SRP_CanonicalAccel(double P0, double earth_radius)
    : P0_(P0), earth_radius_(earth_radius) {}

Eigen::Vector3d SRP_CanonicalAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    Epoch t
) const {
    double current_time = t.toJD();
    Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);

    // Vector from satellite to sun
    Eigen::Vector3d sat_to_sun = sun_pos - pos;
    double dist_sat_to_sun = sat_to_sun.norm();
    Eigen::Vector3d sun_dir = sat_to_sun.normalized();

    // Distance Scaling: Scale flux based on inverse square law
    double au_ratio = AU_ / dist_sat_to_sun;
    double scaled_P0 = P0_ * (au_ratio * au_ratio);

    // Conical Shadow (Penumbra) Math
    double r_sat = pos.norm();
    
    // Apparent angular radii of the Sun and Earth from the satellite
    double a_sun = std::asin(R_SUN_ / dist_sat_to_sun);
    double a_earth = std::asin(earth_radius_ / r_sat);

    // Angle between Earth center and Sun center
    Eigen::Vector3d to_earth = -pos.normalized();
    double theta = std::acos(to_earth.dot(sun_dir));

    double sh = 1.0; // Illumination fraction (1.0 = full sun)

    // Check eclipse conditions
    if (theta <= a_earth - a_sun) {
        // Total eclipse (Umbra)
        sh = 0.0;
    }
    else if (theta < a_earth + a_sun) {
        // Partial eclipse (Penumbra) - Overlapping circles area calculation
        double c = (a_earth * a_earth - a_sun * a_sun + theta * theta) / (2.0 * theta);

        // Protect against domain errors in acos
        double arg1 = std::clamp(c / a_earth, -1.0, 1.0);
        double arg2 = std::clamp((theta - c) / a_sun, -1.0, 1.0);

        double area = a_earth * a_earth * std::acos(arg1) + 
                      a_sun * a_sun * std::acos(arg2) - 
                      a_earth * c * std::sin(std::acos(arg1));

        // Area of the sun disk
        double sun_area = M_PI * a_sun * a_sun;
        sh = 1.0 - (area / sun_area);

        // Clamp to physical bounds just in case of float precision issues
        sh = std::clamp(sh, 0.0, 1.0);
    }

    // Final Acceleration Vector
    return (scaled_P0 * sc.Cr() * sc.area() / sc.mass()) * sh * sun_dir;
}

// =============================
// RelativityAccel
// =============================

RelativityAccel::RelativityAccel(double mu, double c) 
    : mu_(mu), c_(c) {}

Eigen::Vector3d RelativityAccel::computeAcceleration(
    const Spacecraft& sc, 
    const Eigen::Vector3d& pos, 
    const Eigen::Vector3d& vel, 
    Epoch t) const 
{
    // Calculate magnitudes and dot product
    double pos_mag = pos.norm();
    double vel_mag_sq = vel.squaredNorm(); // v^2
    double pos_dot_vel = pos.dot(vel);
    
    double c_sq = c_ * c_;
    
    // Calculate the scalar term scaling the position vector: (4*mu/r - v^2)
    double term1 = (4.0 * mu_ / pos_mag) - vel_mag_sq;
    
    // Full IERS equation// accel = (mu / (c^2 * r^3)) * [ term1 * r + 4 * (r dot v) * v ]
    Eigen::Vector3d accel = (mu_ / (c_sq * pos_mag * pos_mag * pos_mag)) * (term1 * pos + 4.0 * pos_dot_vel * vel);
                            
    return accel;
}

// =============================
// SimpleAlbedoAccel
// =============================

SimpleAlbedoAccel::SimpleAlbedoAccel(double earth_radius)
    : earth_radius_(earth_radius) {}

Eigen::Vector3d SimpleAlbedoAccel::computeAcceleration(
    const Spacecraft& sc,
    const Eigen::Vector3d& pos,
    const Eigen::Vector3d&,
    Epoch t
) const {
    double current_time = t.toJD();
    Eigen::Vector3d sun_pos = SunEphemeris::getSunPositionECI(current_time);

    double cos_sun = sun_pos.dot(pos) / (sun_pos.norm() * pos.norm());

    if (cos_sun <= 0) return Eigen::Vector3d::Zero();

    // F = view factor for a diffusely reflecting sphere
    double F = (earth_radius_ * earth_radius_) / (pos.norm() * pos.norm());
    
    // Albedo Flux in W/m^2
    double albedo_flux = solar_const_ * albedo_coeff_ * F * cos_sun;

    double rad_press = albedo_flux / c_;

    return pos.normalized() * rad_press * sc.Cr() * sc.area() / sc.mass();
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
