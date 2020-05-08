//Disclaimer:
//Some of the function implementations were made with the help of the 
//"Numerical Recipes 3rd Edition: The Art of Scientific Computing"
//by William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery

#include "Numerics.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <limits>       // std::numeric_limits

//1. Simple functions
	int Sign(double arg)
	{
		if(arg > 0.0) 		return 1;
		else if(arg == 0.0)	return 0;
		else 				return -1;
	}

	double Sign(double x, double y)
	{
		if(Sign(x) == Sign(y)) return x;
		else return -1.0 * x;
	}

	double StepFunction(double x)
	{
		if(x>=0) 	return 1.0;
		else		return 0.0;
	}

	double Round(double N,unsigned int digits)
	{
		if(N == 0) return 0;
		unsigned int digits_max = 7;
		if(digits > digits_max)
		{
			std::cerr <<"Error in Round(): Significant digits > "<<digits_max <<"."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		//Make the argument a positive number.
		double sign = Sign(N);
		N *= sign;
		
		//Cut off the decimal power
		double DecimalPower = floor( log10(N) );
		double prefactor = N*pow(10,-DecimalPower);

		//Round the prefactor
		prefactor = std::floor(prefactor * pow(10.0,digits-1)+0.5);
		prefactor = prefactor * pow(10.0,-1.0*digits+1);

		return sign*prefactor*pow(10,DecimalPower);
	}

	double Relative_Difference(double a,double b)
	{
		double d = std::fabs(a-b);
		double max = std::max(fabs(a),fabs(b));
		return d / max;
	}

	bool Floats_Equal(double a, double b, double tol)
	{
		if( Relative_Difference(a,b) < tol) return true;
		else return false;
	}
//2. Special functions
//2.1 Gamma functions
	std::vector<double> FactorialList= {1.0};
	double Factorial(unsigned int n)
	{
		if (n > 170)
		{

			std::cerr <<"Error in Factorial: Overflow for " <<n <<"!."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(n < FactorialList.size()) return FactorialList[n];
		else
		{
			while(FactorialList.size() <= n)
			{
				FactorialList.push_back( FactorialList.back()*FactorialList.size());
			}
			return FactorialList.back();
		}
	}

	double Binomial_Coefficient(int n,int k)
	{
		if(k < 0 || n < 0)
		{
			std::cerr <<"Warning in BinomialCoefficient(): negative arguments. Return 0." <<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(n < k)	return 0;
		else if (n > 170)
		{
			return floor( 0.5 + exp(GammaLn(n+1.0)-GammaLn(k+1.0)-GammaLn(n-k+1.0)) );
		}
		else return floor( 0.5 + Factorial(n)/Factorial(k)/Factorial(n-k) );
	}

	//Logarithmic gamma function
	double GammaLn(double x)
	{
		double cof[14] = {57.1562356658629235,-59.5979603554754912, 14.1360979747417471,-0.491913816097620199,.339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3, -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3, .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
		if(x <= 0)
		{
			std::cerr<<"Error in GammaLn(x): x<=0."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		double sum = 0.999999999999997092;
		double y = x;
		double tmp = x + 671.0 / 128.0;
		tmp = (x+0.5) * log(tmp) - tmp;
		for(int j = 0; j < 14; j++)
		{
			sum += cof[j] / ++y;
		}
		return tmp + log(2.5066282746310005*sum/x);
	}

	double Gamma(double x)
	{
		return exp(GammaLn(x));
	}

	double Upper_Incomplete_Gamma(double x, double s)
	{
		return Gamma(s) * GammaQ(x,s);
	}

	double Lower_Incomplete_Gamma(double x, double s)
	{
		return Gamma(s) * GammaP(x,s);
	}


	//Q(x,a) via integration;
	double GammaQint(double x,double a)
	{
		//Compute P(x,a)
			double gammaP;
			double gln = GammaLn(a);
			//How far to integrate N sqrt(a) around the peak at a-1:
				double N = 10;
				double tPeak = a-1.0;
				double tMin = std::max(0.0, tPeak - N*sqrt(a));
				double tMax = tPeak + N*sqrt(a);
				//If x lies far away from the peak:
				if(x>tMax) gammaP = 1.0;
				else if (x<tMin ) gammaP = 0.0;
			//Numerical integration
				else
				{
					
					//integrand
						std::function<double(double)> integrand = [a,gln] (double t)
						{
							return exp(-gln-t+log(t)*(a-1.0));
						};
						if(x<tMin) tMin = 0.0;
					//Precision
						double eps = Find_Epsilon(integrand,tMin,x,1e-5);
					//Integrate
						gammaP = Integrate(integrand,tMin,x,eps);
				}
		
		return 1.0-gammaP;
	}

	//Series expansion of P(x,a)
	double GammaPser(double x,double a)
	{
		double eps = std::numeric_limits<double>::epsilon();
		double sum,del,ap;
		double gln = GammaLn(a);
		ap = a;
		del = sum = 1.0/a;
		while(fabs(del)>fabs(sum)*eps)
		{
			ap++;
			del *= x/ap;
			sum += del;
		}
		return sum * exp(-x + a*log(x) - gln);
	}

	//Continued fraction representation of Q(x,a)
	double GammaQcf(double x,double a)
	{
		//Precision
			double eps = std::numeric_limits<double>::epsilon();
			double FPMIN = std::numeric_limits<double>::min()/eps;
		double del = 0.0;
		double gln = GammaLn(a);
		double b = x+1.0-a;
		double c = 1.0/FPMIN;
		double d = 1.0/b;
		double h = d;
		int i = 1;
		while(fabs(del-1.0) > eps)
		{
			double an = -1.0*i*(i-a);
			b += 2.0;
			d = an*d+b;
			if(fabs(d) < FPMIN) d = FPMIN;
			c = b + an / c;
			if(fabs(c)<FPMIN) c = FPMIN;
			d = 1.0/d;
			del = d*c;
			h *= del;
		}
		return exp(-x + a*log(x) - gln) * h;
	}

	//Final function using different methods for different parts of the domain
	double GammaQ(double x,double a)
	{
		double aMax = 100.0;
		if(x < 0.0 || a <= 0.0)
		{
			std::cerr <<"Error in GammaQ("<<x<<","<<a<<"): Invalid arguments."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(x == 0) return 1.0;
		else if (a > aMax) return GammaQint(x,a);
		else if (x < a+1.0) return 1.0 - GammaPser(x,a);
		else return GammaQcf(x,a);
	}

	double GammaP(double x,double a)
	{
		return 1.0 - GammaQ(x,a);
	}
				
	//Inverse incomplete gamma function. (Solves P(x,a)=p for x.)
	double Inv_GammaP(double p,double a)
	{
		//Check the arguments
		if(a <= 0.0)
		{
			std::cerr <<"Error in Inv_GammaP(): a must be positive."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		if(p >= 1.0) return std::max(100.0,a + 100.*sqrt(a));
		if(p <= 0.0) return 0.0;

		//Parameter
		double x;
		double gln = GammaLn(a);
		double a1 = a - 1.0;
		double lna1 = log(a1);
		double afac = exp(a1*(lna1-1.0) - gln);

		//Initial guess 1
		if(a > 1.0)
		{
			double pp = (p<0.5)? p : 1.0-p;
			double t = sqrt(-2.0*log(pp));
			x = (2.30753 + t*0.27061) / (1. + t*(0.99229 + t*0.04481)) - t;
			if(p < 0.5) x = -x;
			x = std::max( 1.0e-3, a * pow(1.0 - 1.0/(9.*a) - x/(3.*sqrt(a)),3.0) );

		}
		//Initial guess 2
			else
			{
				double t = 1.0 - a*(0.253+a*0.12);
		        if (p < t) x = pow(p/t,1./a);
		        else x = 1.-log(1.-(p-t)/(1.-t));
			}
		//Halley's method
			double EPS = 1.0e-8;
			for(int i = 0; i < 12; i++)
			{
				if(x <= 0.0) return 0.0;
				double error = GammaP(x,a) - p;
				double t;
				if(a>1.0) t = afac*exp(-(x-a1)+a1*(log(x)-lna1));
				else t = exp(-x+a1*log(x)-gln);
				double u = error/t;
				x -= (t = u/(1.-0.5*std::min(1.,u*((a-1.)/x - 1))));
				if (x <= 0.) x = 0.5*(x + t);
				if (fabs(t) < EPS*x ) break;
			}
		return x;
	}
	//Solves Q(x,a)=p for x.
	double Inv_GammaQ(double q,double a)
	{
		return Inv_GammaP(1.0-q,a);
	}

	//2.2 Other special functions
	double Inv_Erf(double p)
	{
		// return inverfc(1.-p);
		if(fabs(p-1.0) < 1e-16)
		{
			std::cerr <<"Warning in Inv_erf(double): The argument p = "<<p <<" is very close to 1.0. Return 10."<<std::endl;
			return 10.0;
		}
		else if(fabs(p) >= 1.0)
		{
			std::cerr <<"Error in Inv_erf(): Invalid argument |p| = |"<<p <<"| > 1"<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			std::function<double(double)> fct = [p] (double x)
			{
				return erf(x) - p;
			};
			double xLeft = -10.0;
			double xRight = 10.0;
			return Find_Root(fct, xLeft, xRight, 1.0e-4);
		}
	}	

//3. Integration via Adaptive Simpson Method
//3.1 One-dimensional integration via adaptive Simpson method 
	//Function to return a reasonable precision.
	double Find_Epsilon(std::function<double(double)> func, double a,double b,double precision)
	{
		double c = (a + b)/2;
		double h = b - a;                                                                  
			double fa = func(a);
			double fb = func(b);
			double fc = func(c);                                                           
			double S = (h/6)*(fa + 4*fc + fb);
			double epsilon = precision*S;
			return epsilon;
	}

	double Adaptive_Simpson_Integration(std::function<double(double)> func, double a,double b, double epsilon,double S,double fa,double fb,double fc,int bottom,bool &warning)
	{
		double c = (a+b)/2;
		double h = b-a;
		double d = (a+c)/2;
		double e = (b+c)/2;
		double fd = func(d);
		double fe = func(e);
		double Sleft = (h/12)*(fa + 4*fd + fc);                                                           
		double Sright = (h/12)*(fc + 4*fe + fb);                                                          
		double S2 = Sleft + Sright; 
		if (bottom <= 0 || fabs(S2 - S) <= 15*epsilon)//15 due to error analysis 
		{
			if(bottom<=0&&fabs(S2 - S) > 15*epsilon) warning=true;
			return S2 + (S2 - S)/15;  
		}                                         
		else
		{
			return Adaptive_Simpson_Integration(func, a, c, epsilon/2, Sleft,  fa, fc, fd, bottom-1,warning)+Adaptive_Simpson_Integration(func, c, b, epsilon/2, Sright, fc, fb, fe, bottom-1,warning); 
		}    
	}

	//Recursive functions for one-dimensional integration
	double Integrate(std::function<double(double)> func, double a,double b, double epsilon,int maxRecursionDepth)
	{
		int sign = +1;
		if(a == b) return 0.0;
		else if(a > b)
		{
			double aux = a;
			a = b;
			b=aux;
			sign = -1;
		}
		double c = (a + b)/2;
		double h = b - a;                                                                  
			double fa = func(a);
			double fb = func(b);
			double fc = func(c); 
			double S = (h/6)*(fa + 4*fc + fb);
			bool warning = false;
			double result =   Adaptive_Simpson_Integration(func, a, b, fabs(epsilon), S, fa, fb, fc, maxRecursionDepth,warning); 
			if(warning)   
			{
				std::cout <<"Warning in Integrate(): Numerical integration on the interval ("<<a<<","<<b<<") did not converge to the desired precision." <<std::endl;
				std::cout <<"\tDesired precision: " <<Round(fabs(epsilon)) <<" Result: " <<Round(result)<<std::endl;
			}
		if(std::isnan(result)) std::cout <<"Warning in Integrate(): Result is nan."<<std::endl;
		else if(std::isinf(result)) std::cout <<"Warning in Integrate(): Result is inf."<<std::endl;
		return sign*result;
	}

//4. Interpolation
//4.1 One-dimensional interpolation
	void Interpolation::Compute_Steffen_Coefficients(std::vector<std::vector<double>>& data, std::vector<double> &a,std::vector<double> &b,std::vector<double> &c,std::vector<double> &d)
	{
		unsigned int N = data.size();

		//Compute the Steffen coefficients for the interpolation
		//1. h and s.
		// double h[N-1],s[N-1];
		std::vector<double> h(N-1), s(N-1);
		for(unsigned int i = 0; i < N-1; i++)
		{
			double x_i = data[i][0];
			double x_ip1 = data[i+1][0];
			
			double y_i = data[i][1];
			double y_ip1 = data[i+1][1];
			h[i] = x_ip1-x_i;
			s[i] = (y_ip1-y_i)/h[i];
		}

		//2. p and dy
		// double dy[N],p[N];
		std::vector<double> dy(N), p(N);
		for(unsigned int i = 0; i < N; i++)
		{
			//First point
			if(i == 0)
			{
				p[i] = s[i]*(1.0+h[i]/(h[i]+h[i+1]))-s[i+1]*h[i]/(h[i]+h[i+1]);
				dy[i] = (Sign(p[i])+Sign(s[i]))*std::min(1.0*fabs(s[i]),0.5*fabs(p[i]));
			}
			//Last point
			else if(i == N-1)
			{
				p[i] = s[i-1]*(1.0+h[i-1]/(h[i-1]+h[i-2]))-s[i-2]*h[i-1]/(h[i-1]+h[i-2]);
				dy[i] = (Sign(p[i])+Sign(s[i-1]))*std::min(1.0*fabs(s[i-1]),0.5*fabs(p[i]));
			}
			//Points in the middle
			else
			{
				p[i] = (s[i-1]*h[i]+s[i]*h[i-1])/(h[i-1]+h[i]);
				dy[i] = (Sign(s[i-1])+Sign(s[i]))*std::min(1.0*fabs(p[i])/2.0 ,std::min(1.0*fabs(s[i]) ,1.0*fabs(s[i-1]) ) );
			}
		}

		//3. a,b,c, and d
		for(unsigned int i = 0; i < N-1; i++)
		{
			a.push_back((dy[i]+dy[i+1]-2.0*s[i])/pow(h[i],2.0));
			b.push_back((3.0*s[i]-2.0*dy[i]-dy[i+1])/h[i]);
			c.push_back(dy[i]);
			d.push_back(data[i][1]);
			if(std::isnan(a.back())||std::isnan(b.back())||std::isnan(c.back())||std::isnan(d.back())) 
			{
				std::cout <<"Warning: Steffen coefficients in interpolation are NAN."<<std::endl;
			}
		}
	}
	unsigned int Interpolation::Bisection(double x,int jLeft,int jRight)
	{
		while((jRight-jLeft)>1)
		{
			int jm = (jRight+jLeft) >>1 ;
			if (x >= TabulatedData[jm][0])	jLeft=jm;
			else jRight = jm;
		}
		return jLeft;
	}
	unsigned int Interpolation::Hunt(double x)
	{
		// 1. Hunting phase returns jd,ju which bracket j
		int dj = 1;
		int jd;
		unsigned int ju;
		//Hunt up
		if(x>TabulatedData[jLast][0])
		{
			jd = jLast;
			ju = jd+dj;
			while(x>TabulatedData[ju][0])
			{
				jd = ju;
				ju += dj;
				//Check if we ran off the range:
				if(ju > N_Data-1)
				{
					ju = N_Data-1;
					break;
				}
				else
				{
					dj += dj;
				}
			}
		}
		//Hunt down
		else if (x<TabulatedData[jLast][0])
		{
			ju = jLast;
			jd = ju-dj;
			while(x<TabulatedData[jd][0])
			{
				ju = jd;
				jd -= dj;
				//Check if we ran off the range:
				if(jd < 0)
				{
					jd = 0;
					break;
				}
				else
				{
					dj += dj;
				}
			}
		}
		else
		{
			return jLast;
		}

		// 2. Bisection phase
		if((ju-jd)>1) jd = Bisection(x,jd,ju);
		
		return jd;
	}
	// Find j such that list[j]<x<list[j+1]
	unsigned int Interpolation::Locate(double x)
	{
		if( ((xDomain[0]-x) > 0.0) || ((x-xDomain[1]) > 0.0) )
		{
			printf("\nError in Interpolation::Locate(): x = %e lies outside the domain [%e,%e].\n\n",x,xDomain[0],xDomain[1]);
			std::exit(EXIT_FAILURE);
		}
		else
		{
			//Use Bisection() or the Hunt method, depending of the last calls were correlated.
			unsigned int j = corr ? Hunt(x): Bisection(x,0,N_Data-1);
			//Check if the points are still correlated.
			corr = (fabs((j-jLast))<10);

			jLast = j;
			return j;
		}
	}
//Constructors
	Interpolation::Interpolation()
	{
		//Generate some data to interpolate zero
		std::vector<std::vector<double>> data;
		data.push_back(std::vector<double> {-1.0,0.0});
		data.push_back(std::vector<double> {0.0,0.0});
		data.push_back(std::vector<double> {+1.0,0.0});
		
		//Define members
		preFactor = 1.0;
		TabulatedData = data;
		N_Data = data.size();
		xDomain = {-1.0,1.0};
		jLast = 0;
		corr=false;
		
		//Compute coefficients
		Compute_Steffen_Coefficients(data,a,b,c,d);
	}
	Interpolation::Interpolation(const std::string& filename,double dim1,double dim2)
	{
		//1. Import data.
		std::ifstream inputfile;
		inputfile.open(filename);
		if (inputfile.good())
		{
	        while (!inputfile.eof())
	        {
	        	double x,y;
	            inputfile >>x;
	            inputfile >>y;
	            TabulatedData.push_back(std::vector<double> {x*dim1,y*dim2});
	        }
	        inputfile.close();
   		}
   		else
   		{
        	std::cerr << "Error in Interpolation(" <<filename<<"): File does not exist."<<std::endl;
    		std::exit(EXIT_FAILURE);
    	}
		
		//2. Sort data.
		std::sort(TabulatedData.begin(), TabulatedData.end());
	
		//3. Interpolate data.
		preFactor = 1.0;
		N_Data = TabulatedData.size();
		jLast = 0;
		corr=false;
		
		//3.1 Find the function's domain:
		std::vector <double> x;
		for(unsigned i = 0; i < TabulatedData.size(); i++) x.push_back(TabulatedData[i][0]);
		xDomain.push_back( *min_element(x.begin(),x.end()));
		xDomain.push_back( *max_element(x.begin(),x.end()));
		
		//3.2 Compute the Steffen coefficients for the interpolation
		Compute_Steffen_Coefficients(TabulatedData,a,b,c,d);

	}
	Interpolation::Interpolation(const std::vector<std::vector<double>>& data,double dim1,double dim2)
	{
		preFactor = 1.0;
		for(unsigned int i = 0; i < data.size(); i++)
		{
			std::vector<double> aux ={data[i][0]*dim1,data[i][1]*dim2};
			TabulatedData.push_back(aux);
		}
		N_Data = TabulatedData.size();
		jLast = 0;
		corr=false;
		
		//Find the function's domain:
		std::vector <double> x;
		for(unsigned i = 0; i < TabulatedData.size(); i++) x.push_back(TabulatedData[i][0]);
		xDomain.push_back( *min_element(x.begin(),x.end()));
		xDomain.push_back( *max_element(x.begin(),x.end()));
		
		//Compute the Steffen coefficients for the interpolation
		Compute_Steffen_Coefficients(TabulatedData,a,b,c,d);
	}

	std::vector<std::vector<double>> Interpolation::Return_Data()
	{
		return TabulatedData;
	}

	std::vector<double> Interpolation::Return_Domain()
	{
		return xDomain;
	}

	double Interpolation::Return_Prefactor()
	{
		return preFactor;
	}

	std::vector<std::vector<double>> Interpolation::Return_Coefficients()
	{
		std::vector<std::vector<double>> output;
		for (unsigned i = 0 ; i < a.size() ; i++)
		{
			std::vector<double> aux {a[i],b[i],c[i],d[i]};
			output.push_back(aux);
		}
		return output;
	}

	void Interpolation::Set_Prefactor(double factor)
	{
		preFactor=factor;
	}

	double Interpolation::Interpolate(double x)
	{
		int j = Locate(x);
		double x_j = TabulatedData[j][0];
		return preFactor*(a[j]*pow(x-x_j,3.0)+b[j]*pow(x-x_j,2.0)+c[j]*(x-x_j)+d[j]);
	}

	void Interpolation::Multiply(double factor)
	{
		preFactor*=factor;
	}

    void Interpolation::Save_Function(std::string filename,unsigned int points)
    {
    	std::ofstream f;
    	f.open(filename);
    	double dx = (xDomain[1]-xDomain[0])/(points-1);
    	for(unsigned int i = 0;i < points; i++)
 		{
 			double x = xDomain[0] + i * dx;
 			f <<x <<"\t" <<Interpolate(x)<<std::endl;
 		}
    	f.close();
    }
//5. Root finding
	//Root finding with Ridder's method
	double Find_Root(std::function<double(double)> func,double xLeft, double xRight,double xAccuracy)
	{
		const int Max_Iterations = 50;
		//1. Check if xLeft<xRight, otherwise swap.
		if(xLeft>xRight)
		{
			double temp = xLeft;
			xLeft = xRight;
			xRight = temp;
		}

		//2. Compute functions at boundary
		double fLeft = func(xLeft);
		double fRight = func(xRight);
		
		//3. Check if xLeft and xRight bracket a root or already yield a root. Also check for NaN's.
		if(std::isnan(fLeft)||std::isnan(fRight))
		{
			std::cerr <<"Error in Find_Root(): Function returns nan at the brackets."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(fLeft*fRight>=0.0)
		{
			if(fLeft==0) return xLeft;
			else if(fRight==0) return xRight;
			else
			{
				std::cerr<<"Error in Find_Root(): f(xLeft)*f(xRight) = ("<<func(xLeft)<<")*("<<func(xRight) <<")>0.0"<<std::endl;
				std::exit(EXIT_FAILURE);
			}
		}

		//4. Ridder's method
		else
		{
			double x1 = xLeft;
			double x2 = xRight;
			double f1 = fLeft;
			double f2 = fRight;
			double result = -9.9e99;
			for(int i=0;i<Max_Iterations;i++)
			{
				//Mid point
					double x3 = (x1+x2)/2.0;

					double f3 = func(x3);
				//New point
					double x4 = x3 +(x3-x1) * Sign(f1-f2)*f3/sqrt(f3*f3-f1*f2);
				//Check if we found the root
					if(fabs(x4-result)<xAccuracy) return x4;
				//Prepare next iteration
					result = x4;
					double f4 = func(x4);
					if(f4==0.0) return result;
					//a) x3 and x4 bracket the root
					if(Sign(f3,f4) != f3)
					{
						x1 = x3;
						f1 = f3;
						x2 = x4;
						f2 = f4;
					}
					//b) x1 and x4 bracket the root
					else if(Sign(f1,f4) != f1)
					{
						x2 = x4;
						f2 = f4;

					}
					//c) x2 and x4 bracket the root
					else if(Sign(f2,f4) != f2)
					{
						x1 = x4;
						f1 = f4;
					}
					else
					{
						std::cerr <<"Error in Find_Root(). Ridder's method does not reach the root."<<std::endl;
						std::exit(EXIT_FAILURE);
					}
			}
			std::cout <<"Warning in Find_Root(): Iterations exceed the maximum. Final value f("<<result<<")="<<func(result)<<std::endl;
			return result;
		}
	}

//6. Minimization
//6.1 Multi-dimensional
	std::vector<double> Minimization::minimize(std::vector<double> &starting_point, const double delta, std::function<double(std::vector<double>)> func)
	{
		std::vector<double> deltas(starting_point.size(),delta);
    	return minimize(starting_point,deltas,func);
	}

	std::vector<double> Minimization::minimize(std::vector<double> &starting_point, std::vector<double> &deltas, std::function<double(std::vector<double>)> func)
	{
		int ndim=starting_point.size(); 
		std::vector<std::vector<double>> pp(ndim+1,std::vector<double>(ndim,0.0));
		for (int i=0; i<ndim+1; i++)
		{
			for (int j=0; j<ndim; j++) pp[i][j] = starting_point[j];
			if (i !=0 ) pp[i][i-1] += deltas[i-1];
		}
    	return minimize(pp,func);
	}

	std::vector<double> Minimization::minimize(std::vector<std::vector<double>> &pp, std::function<double(std::vector<double>)> func)
	{
		const int NMAX = 5000;
		const double TINY=1.0e-10;
		int ihi,ilo,inhi;
		mpts=pp.size(); //rows
		ndim=pp[0].size(); //columns
		std::vector<double> psum(ndim),pmin(ndim),x(ndim);
		current_simplex=pp;
		y.resize(mpts);
		for (int i=0; i<mpts; i++)
		{
			for (int j=0; j<ndim; j++)
			{
				x[j]=current_simplex[i][j];
			}
			y[i] = func(x);
			
		}

		nfunc=0;
		get_psum(current_simplex,psum);
		for (;;)
		{

			ilo=0;
			//First we must determine which point is the highest (worst), next-highest, and lowest (best), by looping over the points in the simplex.
			ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
			for (int i=0; i<mpts; i++)
			{
				if (y[i] <= y[ilo]) ilo=i;
				if (y[i] > y[ihi])
				{
					inhi=ihi;
					ihi=i;
				}
				else if (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			double rtol=2.0 * fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);
			//Compute the fractional range from highest to lowest and return if satisfactory.
			if (rtol < ftol)
			{
				std::swap(y[0],y[ilo]);
				for (int i=0; i<ndim; i++) 
				{
					std::swap(current_simplex[0][i],current_simplex[ilo][i]);
					pmin[i]=current_simplex[0][i];
				}
				fmin=y[0];
			    return pmin;
			}
			if (nfunc >= NMAX) 
			{
				std::cerr <<"Error in Minimization::minimize(): NMAX exceeded."<<std::endl;
				std::exit(EXIT_FAILURE);
			}
			nfunc += 2;
			
			// Begin a new iteration. First extrapolate by a factor 􏰱1 through the face of the simplex across from the high point, i.e., reflect the simplex from the high point.
			double ytry=amotry(current_simplex,y,psum,ihi,-1.0,func);
			if (ytry <= y[ilo])
			{
				ytry=amotry(current_simplex,y,psum,ihi,2.0,func); //Gives a result better than the best point, so try an additional extrapolation by a factor 2.
			}
			else if (ytry >= y[inhi])
			{
				//The reflected point is worse than the second-highest, so look for an interme- diate lower point, i.e., do a one-dimensional contraction.
				double ysave=y[ihi];
				ytry=amotry(current_simplex,y,psum,ihi,0.5,func);
				if (ytry >= ysave)
				{
					//Can’t seem to get rid of that high point.
					for (int i=0; i<mpts; i++)
					{
						if (i != ilo) 
						{
							for (int j=0; j<ndim; j++) current_simplex[i][j]=psum[j]=0.5*(current_simplex[i][j]+current_simplex[ilo][j]);
							y[i]=func(psum);
						}
					}
					nfunc += ndim; //Keep track of function evaluations.
					get_psum(current_simplex,psum); //Recompute psum.
				}
			}
			else --nfunc; //Correct the evaluation count.
		} // Go back to the test of doneness and the next iteration
	}

	void Minimization::get_psum(std::vector<std::vector<double>> &p, std::vector<double> &psum) //Utility function.
	{
		for (int j=0; j<ndim; j++)
		{
			double sum = 0.0;
			for (int i=0; i<mpts; i++)
			{
				sum += p[i][j];
			}
			psum[j]=sum;
		}
	}

	double Minimization::amotry(std::vector<std::vector<double>> &p, std::vector<double> &y, std::vector<double> &psum,const int ihi, const double fac, std::function<double(std::vector<double>)> func)
	{
		std::vector<double> ptry(ndim);
		double fac1=(1.0-fac)/ndim;
		double fac2=fac1-fac;
		for(int j=0; j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
		double ytry=func(ptry); //Evaluate the function at the trial point.
		if (ytry < y[ihi]) //If it’s better than the highest, then replace the highest.
		{
			y[ihi] = ytry;
			for(int j=0; j<ndim; j++)
			{
				psum[j] += ptry[j]-p[ihi][j];
				p[ihi][j] = ptry[j];
			}
		}
        return ytry;
	}

