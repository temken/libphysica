#ifndef __Natural_Units_hpp_
#define __Natural_Units_hpp_

#include "Linear_Algebra.hpp"

// 1. SI-prefixes
	extern const double yotta, zetta, exa, peta, tera, giga, mega, kilo, hecto, deca, deci, centi, milli, micro, nano, pico, femto, atto, zepto, yocto;

//2. Units
	//2.1 Angles
	extern const double deg, arcmin, arcsec;

	//2.2 Energy
	extern const double GeV, meV, eV, keV, MeV, TeV, PeV, Joule, erg, Rydberg, cal;

	//2.3 Mass
	extern const double gram, kg, tonne, lbs, AMU;

	//2.4 Length
	extern const double cm, mm, meter, km, fm, inch, foot, yard, mile, Angstrom, Bohr_Radius;

	//2.5 Area
	extern const double barn, pb, acre, hectare;

	//2.6 Time
	extern const double sec, ms, ns, minute, hr, day, week, yr;

	//2.7 Frequency
	extern const double Hz;

	//2.8 Force
	extern const double Newton, dyne;

	//2.9 Power
	extern const double Watt;

	//2.10 Pressure
	extern const double Pa, hPa, kPa, bar, barye;

	//2.11 Temperature
	extern const double Kelvin;

	//2.12 Electric Charge
	extern const double Elementary_Charge, Coulomb;

	//2.13 Voltage
	extern const double Volt;

	//2.14 Electric current
	extern const double Ampere;

	//2.15 Electrical capacitance
	extern const double Farad;

	//2.16 Magnetic induction/ flux density
	extern const double Tesla, Gauss;

	//2.17 Magnetic flux
	extern const double Weber;

	//2.20 Electrical resistance
	extern const double Ohm, Siemens;

	//2.21 Amount
	extern const double mole;

//3. Physical constants
	//3.1 Hadron masses
	extern const double mProton, mNeutron, mNucleon;

	//3.2 Quark masses
	extern const double mUp, mDown, mCharm, mStrange, mTop, mBottom;

	//3.3 Lepton masses
	extern const double mElectron, mMuon, mTau;

	//3.4 Boson masses
	extern const double mZ, mW, mHiggs;

	//3.5 Coupling constants
	extern const double aEM, mPlanck, mPlanck_reduced, G_Newton, G_Fermi;
	//3.6 Energy scales
	extern const double Higgs_VeV, QCD_scale;

//4. Astronomical parameters
	//4.1 Masses
	extern const double mEarth;
	extern const double mSun;

	//4.2 Distances
	extern const double rEarth, rSun, AU, pc, kpc, Mpc, ly;

//5. Transform quantities to a given dimension
	extern double In_Units(double quantity, double dimension);
	extern std::vector<double> In_Units(const std::vector<double>&  quantities, double dimension);
	extern std::vector<std::vector<double>> In_Units(const std::vector<std::vector<double>>&  quantities, double dimension);
	extern std::vector<std::vector<double>> In_Units(const std::vector<std::vector<double>>&  quantities, std::vector<double> dimensions);
	extern Vector In_Units(const Vector& quantities, double dimension);
	extern Matrix In_Units(const Matrix& quantities, double dimension);

//6. Simple physics functions
	extern double Reduced_Mass(double m1, double m2);

#endif