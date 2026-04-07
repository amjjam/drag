#ifndef _SUN_MOON_EPHEMERIS_H_
#define _SUN_MOON_EPHEMERIS_H_

#include <cmath>
#include <ctime>
#include <Eigen/Dense>

class SunEphemeris {
public:
    // Calculates Sun position in ECI [meters] given a UTC timestamp
    static Eigen::Vector3d getSunPositionECI(double jd) {

	const double deg2rad = M_PI / 180.0;
	Eigen::Vector3d r_Sun;

        // Days since J2000.0
        const double n = jd - 2451545.0;

        // Mean longitude (L) and anomaly (g) in degrees
        double L = std::fmod(280.460 + 0.9856474 * n, 360.0);
        double g = std::fmod(357.528 + 0.9856003 * n, 360.0);
        if (L < 0) L += 360.0;
        if (g < 0) g += 360.0;

        const double L_rad = L * deg2rad;
        const double g_rad = g * deg2rad;

        // Ecliptic longitude (lambda)
        const double lambda = L + 1.915 * std::sin(g_rad) + 0.020 * std::sin(2 * g_rad);
        const double lambda_rad = lambda * deg2rad;

        // Obliquity of the ecliptic (epsilon)
        const double epsilon_rad = (23.439 - 0.0000004 * n) * deg2rad;

        // Convert to Cartesian ECI coordinates
	// Distance Sun-Earth in AU
	const double r_AU = 1.00014 - 0.01671 * cos(g_rad) - 0.00014 * cos(2 * g_rad);
        // 1 AU = 149,597,870,700 meters
        const double r = r_AU * 149597870700.0;

        r_Sun = r * Eigen::Vector3d(
            			std::cos(lambda_rad),
            			std::cos(epsilon_rad) * std::sin(lambda_rad),
            			std::sin(epsilon_rad) * std::sin(lambda_rad)
        );

	return r_Sun;
    }
};

// Local funtions

inline double Frac(double x) {return x - std::floor(x); }

class MoonEphemeris {
public:
    // Calculates Moon position in ECI [meters] given a Mjd_TT timestamp
    static Eigen::Vector3d getMoonPositionECI(double Mjd_TT) {

	const double deg2rad = M_PI / 180.0;
	const double epsilon_rad = 23.43929111 * deg2rad;
	const double MJD_J2000 = 51544.5;             		// Modif. Julian Date of J2000
	const double Arcs = 3600.0 * 180.0 / M_PI;		// Arcseconds per radian

	const double T   = (Mjd_TT - MJD_J2000) / 36525.0;	// Julian cent. since J2000

  	double  L_0, l,lp, F, D, dL, S, h, N;
  	double  L, B, r, cosB;
  	Eigen::Vector3d  r_Moon;


  	// Mean elements of lunar orbit

  	L_0 =     Frac ( 0.606433 + 1336.851344*T );     // Mean longitude [rev]
                                                   	 // w.r.t. J2000 equinox
  	l   = 2 * M_PI * Frac ( 0.374897 + 1325.552410*T );     // Moon's mean anomaly [rad]
  	lp  = 2 * M_PI * Frac ( 0.993133 +   99.997361*T );     // Sun's mean anomaly [rad]
  	D   = 2 * M_PI * Frac ( 0.827361 + 1236.853086*T );     // Diff. long. Moon-Sun [rad]
  	F   = 2 * M_PI * Frac ( 0.259086 + 1342.227825*T );     // Argument of latitude


  	// Ecliptic longitude (w.r.t. equinox of J2000)

  	dL = +22640*sin(l) - 4586*sin(l-2*D) + 2370*sin(2*D) +  769*sin(2*l)
       	     -668*sin(lp) - 412*sin(2*F) - 212*sin(2*l-2*D) - 206*sin(l+lp-2*D)
       	     +192*sin(l+2*D) - 165*sin(lp-2*D) - 125*sin(D) - 110*sin(l+lp)
       	     +148*sin(l-lp) - 55*sin(2*F-2*D);

  	L = 2 * M_PI * Frac( L_0 + dL/1296.0e3 );  // [rad]

  	// Ecliptic latitude

  	S  = F + (dL+412*sin(2*F)+541*sin(lp)) / Arcs;
  	h  = F-2*D;
  	N  = -526*sin(h) + 44*sin(l+h) - 31*sin(-l+h) - 23*sin(lp+h)
       	     +11*sin(-lp+h) - 25*sin(-2*l+F) + 21*sin(-l+F);

  	B = ( 18520.0*sin(S) + N ) / Arcs;   // Ecliptic Latitude [rad]

  	cosB = cos(B);

  	// Distance [m]

  	r = 385000e3 - 20905e3*cos(l) - 3699e3*cos(2*D-l) - 2956e3*cos(2*D)
     	    -570e3*cos(2*l) + 246e3*cos(2*l-2*D) - 205e3*cos(lp-2*D)
      	    -171e3*cos(l+2*D) - 152e3*cos(l+lp-2*D);

  	// Equatorial coordinates

	Eigen::Matrix3d R_x;
	R_x <<  1,  0,	       0,
		0,  cos(-epsilon_rad), sin(-epsilon_rad),
		0, -sin(-epsilon_rad), cos(-epsilon_rad);

  	r_Moon = r * R_x * Eigen::Vector3d ( cos(L)*cosB, sin(L)*cosB, sin(B) );

  	return r_Moon;

    }
};

#endif
