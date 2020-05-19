#include "Natural_Units.hpp"

#include <cmath>

namespace libphysica
{
namespace natural_units
{

// 1. SI-prefixes
	const double yotta	=	1.0e24;
	const double zetta	=	1.0e21;
	const double exa	=	1.0e18;
	const double peta	=	1.0e15;
	const double tera	=	1.0e12;
	const double giga	=	1.0e9;
	const double mega	=	1.0e6;
	const double kilo	=	1.0e3;
	const double hecto	=	1.0e2;
	const double deca	=	1.0e1;
	const double deci	=	1.0e-1;
	const double centi	=	1.0e-2;
	const double milli	=	1.0e-3;
	const double micro	=	1.0e-6;
	const double nano	=	1.0e-9;
	const double pico	=	1.0e-12;
	const double femto	=	1.0e-15;
	const double atto	=	1.0e-18;
	const double zepto	=	1.0e-21;
	const double yocto	=	1.0e-24;

//2. Units
	//2.1 Angles
	const double deg	=	M_PI / 180.0;
	const double arcmin	=	deg / 60.0;
	const double arcsec	=	arcmin / 60.0;

	//2.2 Energy
	const double GeV		=	1.0;
	const double meV		=	1.0e-12 * GeV;
	const double eV			=	1.0e-9 * GeV;
	const double keV		=	1.0e-6 * GeV;
	const double MeV		=	1.0e-3 * GeV;
	const double TeV		=	1.0e3 * GeV;
	const double PeV		=	1.0e6 * GeV;
	// const double Joule		=	1.0 / (1.602176634e-10) * GeV;
	const double Joule		=	kg*pow(meter/sec,2);
	const double erg		=	gram*pow(cm/sec,2);
	const double Rydberg	=	13.605693009 * eV;
	const double cal		=	4.184 * Joule;

	//2.3 Mass
	const double gram	=	5.60958884493318e23 * GeV;
	const double kg		=	1.0e3 * gram;
	const double tonne	=	1.0e3 * kg;
	const double lbs	=	0.453592 * kg;
	const double AMU	=	0.9314940954 * GeV;

	//2.4 Length
	const double cm				=	5.067730214314311e13 / GeV;
	const double mm				=	1.0e-1 * cm;
	const double meter			=	1.0e2 * cm;
	const double km				=	1.0e3 * meter;
	const double fm				=	1.0e-15 * meter;
	const double inch			=	2.54 * cm;
	const double foot			=	12.0 * inch;
	const double yard			=	3.0 * foot;
	const double mile			=	1609.344 * meter;
	const double Angstrom		=	1.0e-10 * meter;
	const double Bohr_Radius	=	5.291772083e-11 * meter;

	//2.5 Area
	const double barn		=	1.0e-24 * cm*cm;
	const double pb			=	pico * barn;
	const double acre		=	4046.86 * meter*meter;
	const double hectare	=	1.0e4 * meter*meter;

	//2.6 Time
	const double sec = 299792458.0 * meter;
	const double ms = milli * sec;
	const double ns = nano * sec;
	const double minute = 60.0 * sec;
	const double hr = 60 * minute;
	const double day = 24 * hr;
	const double week = 7 * day;
	const double yr = 365.25 * day;

	//2.7 Frequency
	const double Hz = 1.0 / sec;

	//2.8 Force
	const double Newton = kg * meter / sec / sec;
	const double dyne = 1.0e-5 * Newton;

	//2.9 Power
	const double Watt = Joule / sec;

	//2.10 Pressure
	const double Pa		=	Newton / meter / meter;
	const double hPa	=	hecto * Pa;
	const double kPa	=	kilo * Pa;
	const double bar	=	1.0e5 * Pa;
	const double barye	=	dyne / cm / cm;

	//2.11 Temperature
	const double Kelvin	=	8.6173303e-14 * GeV;

	//2.12 Electric Charge
	const double Elementary_Charge	=	0.30282212;
	const double Coulomb			=	Elementary_Charge / (1.602176565e-19);

	//2.13 Voltage
	const double Volt = Joule / Coulomb;

	//2.14 Electric current
	const double Ampere = Coulomb / sec;

	//2.15 Electrical capacitance
	const double Farad = Coulomb / Volt;

	//2.16 Magnetic induction/ flux density
	const double Tesla	=	(Newton * sec) / (Coulomb * meter);
	const double Gauss	=	1.0e-4 * Tesla;

	//2.17 Magnetic flux
	const double Weber = Tesla * meter * meter ;

	//2.20 Electrical resistance
	const double Ohm		=	Volt / Ampere;
	const double Siemens	=	1.0 / Ohm;

	//2.21 Amount
	const double mole = 6.02214076e23;

//3. Physical constants
	//3.1 Hadron masses
	const double mProton	=	938.2720813 * MeV;
	const double mNeutron	=	939.5654133 * MeV;
	const double mNucleon	=	0.932 * GeV;

	//3.2 Quark masses
	const double mUp		=	2.3 * MeV;
	const double mDown		=	4.8 * MeV;
	const double mCharm		=	1.275 * GeV;
	const double mStrange	=	95 * MeV;
	const double mTop		=	173.210 * GeV;
	const double mBottom	=	4.180 * GeV;

	//3.3 Lepton masses
	const double mElectron	=	0.5109989461 * MeV;
	const double mMuon		=	105.6583745 * MeV;
	const double mTau		=	1776.86 * MeV;

	//3.4 Boson masses
	const double mZ		=	91.1876 * GeV;
	const double mW		=	80.379 * GeV;
	const double mHiggs	=	125.18 * GeV;

	//3.5 Coupling constants
	const double aEM				=	1.0/137.035999139;
	const double mPlanck			=	1.22091e19 * GeV;
	const double mPlanck_reduced	=	mPlanck / sqrt(8.0 * M_PI);
	const double G_Newton			=	1.0 / mPlanck / mPlanck;
	const double G_Fermi			=	1.1663787e-5 / GeV / GeV;
	
	//3.6 Energy scales
	const double Higgs_VeV	=	pow((sqrt(2) * G_Fermi),-0.5);
	const double QCD_scale	=	218 * MeV;

//4. Astronomical parameters
	//4.1 Masses
	const double mEarth	=	5.9724e24 * kg;
	const double mSun	=	1.98848e30 * kg;

	//4.2 Distances
	const double rEarth	=	6371 * km;
	const double rSun	=	6.957E8 * meter;
	const double AU		=	149597870700 * meter; 
	const double pc		=	3.08567758e16 * meter;
	const double kpc	=	kilo * pc;
	const double Mpc	=	mega * pc;
	const double ly		=	365.25 * day;

//5. Functions
	double In_Units(double quantity, double dimension)
	{
		return quantity/dimension;
	}

	std::vector<double> In_Units(const std::vector<double>&  quantities, double dimension)
	{
		std::vector<double> result(quantities.size());
		for(unsigned int i = 0; i < quantities.size(); i++)
		{
			result[i] = In_Units(quantities[i],dimension);
		}
		return result;
	}
	std::vector<std::vector<double>> In_Units(const std::vector<std::vector<double>>&  quantities, double dimension)
	{
		std::vector<std::vector<double>> result(quantities.size());
		for(unsigned int i = 0; i < quantities.size(); i++)
		{
			result[i] = In_Units(quantities[i],dimension);
		}
		return result;
	}

	std::vector<std::vector<double>> In_Units(const std::vector<std::vector<double>>&  quantities, std::vector<double> dimensions)
	{
		std::vector<std::vector<double>> result(quantities.size(), std::vector<double>(dimensions.size(),0.0));
		for(unsigned int i = 0; i < quantities.size(); i++)
		{
			if(quantities[i].size() != dimensions.size())
			{
				std::cerr <<"Error in In_Units(const std::vector<std::vector<double>>&,std::vector<double>): Dimensions of arrays are not consistent."<<std::endl;
				std::exit(EXIT_FAILURE);
			}
			else
			{
				for(unsigned int j = 0; j < quantities[i].size(); j++)
				{
					result[i][j] = In_Units(quantities[i][j],dimensions[j]);
				}
			}
		}
		return result;
	}

	Vector In_Units(const Vector& quantities, double dimension)
	{
		Vector result(quantities.Size());
		for(unsigned int i = 0; i < quantities.Size(); i++)
		{
			result[i] = In_Units(quantities[i],dimension);
		}
		return result;
	}


	Matrix In_Units(const Matrix& quantities, double dimension)
	{
		return 1.0/dimension * quantities;
	}

}	// namespace natural_units

//6. Simple physics functions
	double Reduced_Mass(double m1,double m2)
	{
		return m1*m2/(m1+m2);
	}
	
}	// namespace libphysica
