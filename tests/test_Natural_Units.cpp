#include "gtest/gtest.h"

#include <cmath>

#include "libphysica/Natural_Units.hpp"

using namespace libphysica;
using namespace libphysica::natural_units;

//2. Units
TEST(TestNaturalUnits, TestUnitsAngles)
{
	ASSERT_DOUBLE_EQ(deg, M_PI / 180.0);
	ASSERT_DOUBLE_EQ(arcmin, deg / 60.0);
	ASSERT_DOUBLE_EQ(arcsec, deg / 60.0 / 60.0);
}

TEST(TestNaturalUnits, TestUnitsEnergy)
{
	ASSERT_DOUBLE_EQ(eV, 1.0e-9 * GeV);
	ASSERT_DOUBLE_EQ(meV, 1.0e-3 * eV);
	ASSERT_DOUBLE_EQ(keV, 1.0e3 * eV);
	ASSERT_DOUBLE_EQ(MeV, 1.0e6 * eV);
	ASSERT_DOUBLE_EQ(TeV, 1.0e12 * eV);
	ASSERT_DOUBLE_EQ(PeV, 1.0e15 * eV);
	ASSERT_DOUBLE_EQ(Joule, kg * meter * meter / sec / sec);
	ASSERT_DOUBLE_EQ(erg, gram * cm * cm / sec / sec);
	ASSERT_DOUBLE_EQ(cal, 4.184 * Joule);
}

TEST(TestNaturalUnits, TestUnitsMass)
{
	ASSERT_DOUBLE_EQ(kg, 1.0e3 * gram);
	ASSERT_DOUBLE_EQ(tonne, 1.0e3 * kg);
}

TEST(TestNaturalUnits, TestUnitsLength)
{
	ASSERT_DOUBLE_EQ(mm, 1.0e-1 * cm);
	ASSERT_DOUBLE_EQ(meter, 1.0e2 * cm);
	ASSERT_DOUBLE_EQ(km, 1.0e5 * cm);
	ASSERT_DOUBLE_EQ(fm, 1.0e-13 * cm);
	ASSERT_DOUBLE_EQ(foot, 12.0 * inch);
	ASSERT_DOUBLE_EQ(yard, 3.0 * foot);
	ASSERT_DOUBLE_EQ(mile, 1609.344 * meter);
}

TEST(TestNaturalUnits, TestUnitsTime)
{
	ASSERT_DOUBLE_EQ(minute, 60.0 * sec);
	ASSERT_DOUBLE_EQ(hr, 60.0 * minute);
	ASSERT_DOUBLE_EQ(day, 60.0 * 24.0 * minute);
	ASSERT_DOUBLE_EQ(year, 365.25 * 60.0 * 24.0 * minute);
}

TEST(TestNaturalUnits, TestUnitsForce)
{
	ASSERT_DOUBLE_EQ(Newton, kg * meter / sec / sec);
	ASSERT_DOUBLE_EQ(dyne, gram * cm / sec / sec);
}

TEST(TestNaturalUnits, TestUnitsPower)
{
	ASSERT_DOUBLE_EQ(Watt, Joule / sec);
}

//5. Transform quantities to a given dimension
TEST(TestNaturalUnits, TestInUnitsScalar)
{
	// ARRANGE
	double c				  = 1.0;
	double speed_of_light_mps = 299792458.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(In_Units(c, meter / sec), speed_of_light_mps);
}

TEST(TestNaturalUnits, TestInUnitsStdVector)
{
	// ARRANGE
	std::vector<double> quantities	   = {3.0 * Ohm, 4.0 * Ohm, 2.5 * Ohm};
	std::vector<double> quantities_Ohm = {3.0, 4.0, 2.5};
	// ACT & ASSERT
	ASSERT_EQ(In_Units(quantities, Ohm).size(), quantities_Ohm.size());
	for(unsigned int i = 0; i < quantities_Ohm.size(); i++)
		ASSERT_DOUBLE_EQ(In_Units(quantities, Ohm)[i], quantities_Ohm[i]);
}

TEST(TestNaturalUnits, TestInUnitsNestedStdVector1)
{
	// ARRANGE
	std::vector<std::vector<double>> quantities	   = {{3.0 * Hz, 4.0 * Hz}, {4.0 * Hz, 5.0 * Hz}};
	std::vector<std::vector<double>> quantities_Hz = {{3.0, 4.0}, {4.0, 5.0}};
	// ACT & ASSERT
	ASSERT_EQ(In_Units(quantities, Hz).size(), quantities_Hz.size());
	for(unsigned int i = 0; i < quantities_Hz.size(); i++)
	{
		ASSERT_EQ(In_Units(quantities, Hz)[i].size(), quantities_Hz[i].size());
		for(unsigned int j = 0; j < quantities_Hz[i].size(); j++)
			ASSERT_DOUBLE_EQ(In_Units(quantities, Hz)[i][j], quantities_Hz[i][j]);
	}
}

TEST(TestNaturalUnits, TestInUnitsNestedStdVector2)
{
	// ARRANGE
	std::vector<std::vector<double>> quantities			= {{0.1 * hr, 1.5 * Watt}, {0.3 * hr, 2.3 * Watt}};
	std::vector<std::vector<double>> quantities_hr_Watt = {{0.1, 1.5}, {0.3, 2.3}};
	std::vector<double> units							= {hr, Watt};
	// ACT & ASSERT
	ASSERT_EQ(In_Units(quantities, units).size(), quantities_hr_Watt.size());
	for(unsigned int i = 0; i < quantities_hr_Watt.size(); i++)
	{
		ASSERT_EQ(In_Units(quantities, {hr, Watt})[i].size(), quantities_hr_Watt[i].size());
		for(unsigned int j = 0; j < quantities_hr_Watt[i].size(); j++)
			ASSERT_DOUBLE_EQ(In_Units(quantities, {hr, Watt})[i][j], quantities_hr_Watt[i][j]);
	}
}

TEST(TestNaturalUnits, TestInUnitsVector)
{
	// ARRANGE
	Vector vec({1.1 * Newton, 2.2 * Newton, 3.3 * Newton});
	Vector vec_Newton({1.1, 2.2, 3.3});
	// ACT & ASSERT
	ASSERT_EQ(In_Units(vec, Newton).Size(), vec_Newton.Size());
	for(unsigned int i = 0; i < vec_Newton.Size(); i++)
		ASSERT_DOUBLE_EQ(In_Units(vec, Newton)[i], vec_Newton[i]);
}

TEST(TestNaturalUnits, TestInUnitsMatrix)
{
	// ARRANGE
	Matrix M({{13.3 * AU, 12.2 * AU}, {11.1 * AU, 10.0 * AU}, {9.9 * AU, 8.8 * AU}});
	Matrix M_AU({{13.3, 12.2}, {11.1, 10.0}, {9.9, 8.8}});

	// ACT & ASSERT
	ASSERT_EQ(In_Units(M, AU).Rows(), M_AU.Rows());
	ASSERT_EQ(In_Units(M, AU).Columns(), M_AU.Columns());
	for(unsigned int i = 0; i < M.Rows(); i++)
		for(unsigned int j = 0; j < M.Columns(); j++)
			ASSERT_DOUBLE_EQ(In_Units(M, AU)[i][j], M_AU[i][j]);
}
//6. Simple physics functions
TEST(TestNaturalUnits, TestReducedMass)
{
	// ARRANGE
	double m1		  = 1.0;
	double m2		  = 2.0;
	double mu_correct = 2.0 / 3.0;
	// ACT & ASSERT
	ASSERT_DOUBLE_EQ(Reduced_Mass(m1, m2), mu_correct);
}
