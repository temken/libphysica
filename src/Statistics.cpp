#include "Statistics.hpp"

#include <iostream>
#include <algorithm>

//1. Distributions
	//1.1. Uniform distribution
	double PDF_Uniform(double x, double x_min, double x_max)
	{
		return 1.0/(x_max-x_min);
	}

	double CDF_Uniform(double x, double x_min, double x_max)
	{
		return (x - x_min) / (x_max - x_min);
	}

	//1.2 Normal distribution
	double PDF_Gauss(double x, double mu, double sigma)
	{
		return 1.0/sqrt(2.0*M_PI)/sigma*exp(-pow((x-mu)/sigma,2.0)/2.0);
	}

	double CDF_Gauss(double x, double mu, double sigma)
	{
		return 0.5 * (1.0 + erf( (x-mu) / (sqrt(2)*sigma)));
	}

	double Quantile_Gauss(double p, double mu, double sigma)
	{
		return mu + sqrt(2.0) * sigma * Inv_erf(2.0*p-1.0);
	}

	//1.3 Binomial distribution
	double PMF_Binomial(int mu, double p, int trials)
	{
		return Binomial_Coefficient(mu,trials)*pow(p,trials)*pow(1.0-p,(mu-trials));;
	}

	double CDF_Binomial(int mu, double p, int trials)
	{
		double cdf = 0.0;
		for(int trial = 0; trial <= trials; trial++) cdf += PMF_Binomial(mu,p,trial);
		return cdf;
	}

	//1.4 Poission distribution
	double PMF_Poisson(double expected_events, unsigned int events)
	{
		if(expected_events==0&&events==0) return 1.0;
		else if(expected_events==0&&events>0) return 0.0;
		else
		{
			double sum=events*log(expected_events) -expected_events;
			for(unsigned int i=2;i<=events;i++)sum-= log(i);
			return exp(sum);
		}
	}

	double CDF_Poisson(double expectation_value,unsigned int observed_events)
	{
		double gq = GammaQ(expectation_value,observed_events+1);
		if (gq>=0) return gq;
		else return 0.0;
	}

	double Inv_CDF_Poisson(unsigned int observed_events, double cdf)
	{
		if(observed_events==0) return (-1.0)*log(cdf);
		else return Inv_GammaQ(cdf,observed_events+1); 
	}

	//1.5 Chi-square distribution
	double PDF_Chi_Square(double x, double dof)
	{
		if(x < 0)
		{
			std::cerr <<"Error in PDF_Chi_Square(double, double): Negative x value: x = "<<x<<"."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else return 1.0/ pow(2.0,dof/2.0) / Gamma(dof/2.0) * pow(x,dof/2.0-1.0) * exp(-x/2.0);
	}

	double CDF_Chi_Square(double x, double dof)
	{
		if(x < 0)
		{
			std::cerr <<"Error in CDF_Chi_Square(double, double): Negative x value: x = "<<x<<"."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(fabs(dof)<1e-6)  return 1.0;
		else return 1.0/Gamma(dof/2.0) * Lower_Incomplete_Gamma(x/2.0,dof/2.0);
	}

	double PDF_Chi_Bar_Square(double x, std::vector<double> weights)
	{
		if(x < 0)
		{
			std::cerr <<"Error in PDF_Chi_Bar_Square(double, std::vector<double>): Negative x value: x = "<<x<<"."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			double pdf = 0.0;
			for(unsigned int dof = 1; dof < weights.size(); dof++) //start at 1 because of dof = 0
			{
				pdf += weights[dof] * PDF_Chi_Square(x,dof);
			}
			return pdf;
		}
	}

	double CDF_Chi_Bar_Square(double x, std::vector<double> weights)
	{
		if(x < 0)
		{
			std::cerr <<"Error in CDF_Chi_Bar_Square(double, std::vector<double>): Negative x value: x = "<<x<<"."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else if(x == 0) return 0.0;
		else
		{
			double cdf = weights[0];
			for(unsigned int dof = 1; dof < weights.size(); dof++) //start at 1 because of dof = 0
			{
				cdf += weights[dof] * CDF_Chi_Square(x,dof);
			}
			if(cdf > 1.0)
			{
				std::cerr<<"Warning in CDF_Chi_Bar_Square(double,std::vector<double>): 1-CDF = "<<(1.0 - cdf)<<"< 0. Return 1." <<std::endl;
				return 1.0;
			}
			else return cdf;
		}
	}

	//1.6 Exponential distribution
	double PDF_Exponential(double x, double mean)
	{
		return 1.0 / mean * exp(-1.0 / mean * x);
	}

	double CDF_Exponential(double x, double mean)
	{
		return 1.0 - exp(-1.0 / mean * x);
	}

//2. Likelihoods
	double Log_Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background)
	{
		double log_N_obs_factorial = 0.0;
		for(int j=1; j<= N_observed; j++) log_N_obs_factorial += log(j);
		
		return N_observed * log(N_prediction + expected_background) - log_N_obs_factorial - (N_prediction + expected_background);
	}

	double Likelihood_Poisson(double N_prediction, unsigned long int N_observed, double expected_background)
	{
		return exp(Log_Likelihood_Poisson(N_prediction, N_observed, expected_background));
	}

	double Log_Likelihood_Poisson_Binned(std::vector<double> N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned )
	{
		unsigned int N_bins = N_prediction_binned.size();
		if(expected_background_binned.empty()) expected_background_binned = std::vector<double>(N_bins, 0.0);
		if(N_observed_binned.size() != N_bins || expected_background_binned.size() != N_bins)
		{
			std::cerr <<"Error in Log_Likelihood_Poisson_Binned(): Predicted signals, observed events, and/or expected background are not of equal size."<<std::endl;
			std::exit(EXIT_FAILURE);
		}
		else
		{
			double log_likelihood = 0.0;
			for(unsigned int i = 0; i < N_bins; i++)
			{
				log_likelihood += Log_Likelihood_Poisson(N_prediction_binned[i], N_observed_binned[i], expected_background_binned[i]);
			}
			return log_likelihood;
		}
	}


	double Likelihood_Poisson_Binned(std::vector<double> N_prediction_binned, const std::vector<unsigned long int>& N_observed_binned, std::vector<double> expected_background_binned)
	{
		return exp( Log_Likelihood_Poisson_Binned(N_prediction_binned, N_observed_binned, expected_background_binned) );
	}


//3. Sampling random numbers
	// Real random number which is used as a seed for the PRNG:
	std::random_device rd;
	//The PRNG:
	std::mt19937 generator(rd());
	//Propability:
	std::uniform_real_distribution<double> Sample_Xi(0,1);

	//3.1 Sample from specific distribution
	double Sample_Uniform(std::mt19937& PRNG, double x_min, double x_max)
	{
		return x_min + Sample_Xi(PRNG) * (x_max - x_min);
	}

	int Sample_Poisson(std::mt19937& PRNG, double expectation_value) //Algorithm from https://en.wikipedia.org/wiki/Poisson_distribution
	{
		double STEP = 500;
		double lambda_left = expectation_value;
		int k = 0;
		double p = 1.0;
		do
		{
			k++;
			double u = Sample_Xi(PRNG);
			p = p*u;
			while(p < 1.0 && lambda_left > 0.0)
			{
				if(lambda_left > STEP)
				{
					p = p * exp(STEP);
					lambda_left -= STEP;
				}
				else
				{
					p = p * exp(lambda_left);
					lambda_left = 0.0;
				}	
			}
		} while (p > 1);
		return (k-1);
	}

	//3.2 General sampling algorithms
	double Rejection_Sampling(const std::function<double(double)>& PDF,double xMin,double xMax,double yMax,std::mt19937& PRNG)
	{
		bool success = false;
		double x;
		int count = 0;
		while(!success)
		{
			count++;
			if(count%1000 == 0)
			{
				std::cout <<"Warning in Rejection_Sampling(): Very inefficient sampling with N="<<count<<"."<< std::endl;
			}
			x = xMin + Sample_Xi(PRNG) * (xMax-xMin);
			double y = Sample_Xi(PRNG) * yMax;
			double pdf = PDF(x);
			if (pdf < 0.0) 
			{
				std::cerr <<"Error in Rejection_Sampling(): PDF is negative -> f("<<x<<")=" <<pdf<< std::endl;
				std::exit(EXIT_FAILURE);
			}
			else if(pdf>yMax)
			{
				std::cout <<"Warning in Rejection_Sampling(): PDF>yMax, yMax is set to PDF."<< std::endl;
				return Rejection_Sampling(PDF,xMin,xMax,pdf,PRNG);
			}
			else if(y <= pdf) success = true;
		}
		return x;
	}
	double Inverse_Transform_Sampling(const std::function<double(double)>& cdf,double xMin,double xMax,std::mt19937& PRNG)
	{
		double xi = Sample_Xi(PRNG);
		std::function<double(double)> fct =[&cdf,xi](double x)
		{
			return xi-cdf(x);
		};
		return Find_Root(fct,xMin,xMax,1e-10*(xMax-xMin));
	}

//4. Data point with statistical weight
 	DataPoint::DataPoint(double v, double w)
    {
    	value=v;
    	weight=w;
    }
        
    bool operator <(const DataPoint& lhs,const DataPoint& rhs)
    {
    	return (lhs.value<rhs.value);
    }
	
	bool operator >(const DataPoint& lhs,const DataPoint& rhs)
	{
		return (lhs.value>rhs.value);
	}
	
	bool operator ==(const DataPoint& lhs,const DataPoint& rhs)
	{
		return (lhs.value==rhs.value);
	}
	
	std::ostream& operator <<(std::ostream &output,const DataPoint& dp)
	{
		output <<dp.value<<"\t" <<dp.weight;
		return output;
	}

//5. Basic data analysis
	double Arithmetic_Mean(const std::vector<double>& data)
	{
		return 1.0 * std::accumulate(data.begin(),data.end(),0.0) / data.size();
	}

	double Median(std::vector<double>& data)
	{
		if (data.size() % 2 == 0)
		{
		    const auto median_it1 = data.begin() + data.size() / 2 - 1;
		    const auto median_it2 = data.begin() + data.size() / 2;
		    std::nth_element(data.begin(), median_it1 , data.end());
		    const auto e1 = *median_it1;
		    std::nth_element(data.begin(), median_it2 , data.end());
		    const auto e2 = *median_it2;
		    return (e1 + e2) / 2;
		}
		else
		{
		    const auto median_it = data.begin() + data.size() / 2;
		    std::nth_element(data.begin(), median_it , data.end());
		    return *median_it;
		}
	}

	double Variance(const std::vector<double>& data)
	{
		double mean = Arithmetic_Mean(data);
		double variance = 0.0;
		for(unsigned int i = 0; i < data.size(); i++)
		{
			variance += (data[i] - mean) * (data[i] - mean);
		}
		variance = 1.0 * variance / data.size();
		return variance;
	}

	double Standard_Deviation(const std::vector<double>& data)
	{
		return sqrt(Variance(data));
	}

	std::vector<double> Weighted_Average(std::vector<DataPoint>& data)
	{
		double sum =0.0;
		double wsum=0.0;
		long unsigned int N=data.size();
		//1. Average
			for(unsigned int i=0;i<N;i++)
			{
				sum+=data[i].weight*data[i].value;
				wsum+=data[i].weight;
			}
			double Average = sum/wsum;
			double wAverage= wsum/N;
		//2. Standard error (Cochran)
			double sum1=0.0,sum2=0.0,sum3=0.0;
			for(unsigned int i=0;i<data.size();i++)
			{
				sum1+=pow((data[i].weight*data[i].value-Average*wAverage),2.0);
				sum2+=(data[i].weight-wAverage)*(data[i].weight*data[i].value-Average*wAverage);
				sum3+=pow(data[i].weight-wAverage,2.0);
			}
			double SE=N/(N-1.0)/wsum/wsum*(sum1-2.0*Average*sum2+pow(Average,2.0)*sum3);
		//3. Return result
			return std::vector<double> {Average,sqrt(SE)};
	}

//6. Kernel density estimation
	double Gaussian_Kernel(double x)
	{
		return PDF_Gauss(x, 0.0, 1.0);
	}
	Interpolation Perform_KDE(std::vector<DataPoint> data,double xMin,double xMax,double bw)
	{
		unsigned int N_Data = data.size();
		double Weight_Sum=0.0;

		//Count the weights
		for(unsigned int i = 0 ; i<N_Data ; i++) Weight_Sum += data[i].weight;

		//1. Bandwidth selection:
		//If the bandwidth is not set manually, we find it here.
			if(bw==0)
			{
				//1.1 Average.
					double AverageSum=0.0;
					for(unsigned int i = 0 ; i<N_Data ; i++)
					{
						AverageSum+=data[i].weight*data[i].value;
					}
					double Average = AverageSum/Weight_Sum;
				//1.2 Standard deviation
					double Variance=0.0;
					for(unsigned int i = 0 ; i<N_Data ; i++)
					{
						Variance += data[i].weight * pow(data[i].value - Average,2.0) / Weight_Sum;
					}
				//1.3 Bandwidth with rule-of-thumb estimator
					bw = sqrt(Variance)*pow(4.0/3.0/N_Data,0.2) ;
			}
		
		//Sort data:
		std::sort(data.begin(), data.end());
		// Pseudo data
		unsigned int N_PseudoData = N_Data/3.0;
		
		//2. Perform and tabulate the KDE
			int points = 150;
			double dx = (xMax-xMin)/(points-1);
			std::vector<std::vector<double>> Interpol_List;
			for(int j =0;j<points;j++)
			{
				double x = xMin + j*dx;
				double kde=0.0;
				for(unsigned int i = 0;i<N_Data;i++)
				{
					//Vanilla Gauss KDE
						kde+= data[i].weight*Gaussian_Kernel((x-data[i].value)/bw);
					//Cowling and Hall pseudo data method 1
						if(i<N_PseudoData)
						{
							double xPseudo = 4.0*xMin-6.0*data[i].value+4.0*data[2*i].value-data[3*i].value;
							double wPseudo = (data[i].weight+data[2*i].weight+data[3*i].weight)/3.0;
							kde+= wPseudo*Gaussian_Kernel((x-xPseudo)/bw);
						}

					//Local kernel normalization close to the boundary
						// double data_point = data[i].value;
						// if(data_point < xMin + 3.0* bw)
						// {
						// 	std::function<double(double)> integrand = [data_point,bw] (double x)
						// 	{
						// 		return Gaussian_Kernel((x-data_point)/bw) / bw;
						// 	};
						// 	double norm = Integrate(integrand,xMin, xMax,1e-5);
						// 	kde+= data[i].weight*Gaussian_Kernel((x-data[i].value)/bw) / norm;
						// 	std::cout <<x <<bw <<"\t" <<norm <<std::endl;
						// }
						// else kde+= data[i].weight*Gaussian_Kernel((x-data[i].value)/bw);

					//Cowling and Hall pseudo data method 2
						// double xi3 = (i%3 == 0)? data[i/3].value : 2.0/3.0 * data[std::floor(i/3.0)].value + 1.0/3.0 * data[std::ceil(i/3.0)].value;
						// double x2i3 = (i%3 == 0)? data[2*i/3].value : 1.0/3.0 * data[std::floor(i/3.0)].value + 2.0/3.0 * data[std::ceil(i/3.0)].value;
						// double xi = data[i].value;

						// double wi3 = (i%3 == 0)? data[i/3].weight : 2.0/3.0 * data[std::floor(i/3.0)].weight + 1.0/3.0 * data[std::ceil(i/3.0)].weight;
						// double w2i3 = (i%3 == 0)? data[2*i/3].weight : 1.0/3.0 * data[std::floor(i/3.0)].weight + 2.0/3.0 * data[std::ceil(i/3.0)].weight;
						// double wi = data[i].weight;

						// double xPseudo = 20.0/3.0 * xMin - 5.0* xi3 - 4.0 * x2i3 + 10.0/3.0 * xi;
						// double wPseudo = (wi3 + w2i3 + wi)/3.0;

						// kde += wPseudo*Gaussian_Kernel((x-xPseudo)/bw);
					//Reflection
						// double xRefl = 2.0 * xMin - data[i].value;
						// kde+=data[i].weight*Gaussian_Kernel((x-xRefl)/bw);
				}				
				kde/= bw*Weight_Sum;
				Interpol_List.push_back( std::vector<double> {x,kde});
			}

		Interpolation result(Interpol_List);

		//3. Check normalization/ re-normalize.
			double norm = Integrate(result,xMin,xMax,1e-5);
			result.Multiply(1.0/norm);

			return result;
	}
