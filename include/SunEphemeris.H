#ifndef _SUN_EPHEMERIS_H_
#define _SUN_EPHEMERIS_H_

#include <cmath>
#include <ctime>
#include <Eigen/Dense>

class SunEphemeris {
public:
    // Calculates Sun position in ECI [meters] given a UTC timestamp
    static Eigen::Vector3d getSunPositionECI(std::time_t t) {
        // 1. Convert to Julian Date (JD)
        // 2440587.5 is the JD for 1970-01-01 00:00:00 UTC
        double jd = 2440587.5 + (static_cast<double>(t) / 86400.0);

        // 2. Days since J2000.0
        double n = jd - 2451545.0;

        // 3. Mean longitude (L) and anomaly (g) in degrees
        double L = std::fmod(280.460 + 0.9856474 * n, 360.0);
        double g = std::fmod(357.528 + 0.9856003 * n, 360.0);
        if (L < 0) L += 360.0;
        if (g < 0) g += 360.0;

        double deg2rad = M_PI / 180.0;
        double L_rad = L * deg2rad;
        double g_rad = g * deg2rad;

        // 4. Ecliptic longitude (lambda)
        double lambda = L + 1.915 * std::sin(g_rad) + 0.020 * std::sin(2 * g_rad);
        double lambda_rad = lambda * deg2rad;

        // 5. Obliquity of the ecliptic (epsilon)
        double epsilon_rad = (23.439 - 0.0000004 * n) * deg2rad;

        // 6. Convert to Cartesian ECI coordinates
        // 1 AU = 149,597,870,700 meters
        double r = 149597870700.0; 

        return Eigen::Vector3d(
            r * std::cos(lambda_rad),
            r * std::cos(epsilon_rad) * std::sin(lambda_rad),
            r * std::sin(epsilon_rad) * std::sin(lambda_rad)
        );
    }
};

#endif
