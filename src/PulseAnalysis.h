/*************************************************************************
 *
 *  Filename: PulseAnalysis.h
 *
 *  Description:
 *
 *	Author(s):
 *     Michael T. Febbraro
 *
 *
 *  Creation Date: 11/25/2012
 *  Last modified: 3/21/2016
 *
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <deque>
#include <vector>


using namespace std;

class PulseAnalysis {
private :

protected :

	string  version;
	string  revision;

	float				deltaT, tempV;
	float 				PSD, integralS, integralL, sumL, sumS;
	float 				linear[4], pulse_temp[2000], summing[2000], MovingAverage[2000];
	float 				w,x,y,z;
	int 				k,l,m,n;
	int 				return_code	;
	bool				positive, negative;

public :

						PulseAnalysis();
						~PulseAnalysis();
	void				GetVersion();
	int	   			PSD_Integration (const std::vector< float >, int, int, int, int, float*, float*);
    int             PSD_Integration_Afterpulsing(float*, int, int, int, int, int, float, int*, float*, float*, int*);
	int				Baseline_restore (std::vector< float > &, float*, int, int);
	int     		PSD_Zerocross (float*, int, int, int, float*);
	int				Parameters (float*, int, int, float*, float*, float*, float*, float*);
	int				Parameters2 (float*, int, int, float*, float*);
	int				Time_Pickoff (float*, int, int, int, int, int, float*);
	int				Derivative(float*, int, int);
	int				Integral(float*, int);
	int				PeakFinder (float*, int, int, int, int, int*, int*);
	int				OptimizePSD (const std::vector< float >, int, int, int, int, float*, float*);
	int				Half_Integral(float*, int, float*);
	int				Smooth(float*, int, int, int, float);
	int				HPGe(float*, int, float*);
	std::vector< float >	CMA_Filter(std::vector<float>, size_t, float, float);
};

/** Constructor  */
PulseAnalysis::PulseAnalysis()
{
	deltaT = 1;
}

/** Destructor  */
PulseAnalysis::~PulseAnalysis()
{

}

/** ----------------------------------------------------
*	Get version
*		-  Print the version and revision date
*
*  Notes:
*		Keep this information up to date.
*	----------------------------------------------------
*/
void PulseAnalysis::GetVersion()
{
	string  version   = "Pulse analysis 1.0";
	string  revision  = "9/17/2013";

	cout << "Package: " << version << endl;
	cout << "Revision: " << revision << endl;
}

/** ----------------------------------------------------
*	Continous moving average (CMA) filter
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			halfWidth - halfwidth for moving average
*			preloadValue - initial guess for baseline
*			rejectThreshold - threshold
*  Returns:
*			CMAtrace - output CMA trace
*	----------------------------------------------------
*/
std::vector< float > PulseAnalysis::CMA_Filter(std::vector< float > waveform, size_t halfWidth, float preloadValue, float rejectThreshold)
{
	deque<float> movingBaselineFilter;
	std::vector< float > CMAtrace;

	float movingBaselineValue = 0;

	if( preloadValue > -7777 ) {
		for(size_t i = 0; i < halfWidth; i++ ) {
			movingBaselineFilter.push_back( preloadValue );
			movingBaselineValue += preloadValue;
		}
	}

	for(size_t i = 0; i < halfWidth; i++ ) {
		if( preloadValue > -7777 ) {
			if( fabs( waveform[i] - movingBaselineValue/movingBaselineFilter.size() ) >= rejectThreshold ) {
				continue;
			}
		}
		movingBaselineFilter.push_back( waveform[i] );
		movingBaselineValue += waveform[i];
	}

	for(auto value : waveform) {
		if( fabs( value - movingBaselineValue/movingBaselineFilter.size() ) < rejectThreshold  ) {
			if( movingBaselineFilter.size() >= halfWidth * 2 + 1 ) {
				// filter is fully occupied, pop as we move
				movingBaselineValue -= movingBaselineFilter.front();
				movingBaselineFilter.pop_front();
			}

			movingBaselineValue += value;
			movingBaselineFilter.push_back( value );
		}
		CMAtrace.push_back(movingBaselineValue / movingBaselineFilter.size());
	}

	return CMAtrace;
}

/** ----------------------------------------------------
*	Optimize PSD by charge integration method
*		- Calculates discrimination parameters by
*		  integration of user defined regions of an
*		  input pulse over a user defined range for evaluation.
*
*	Inputs:
*			pulse - input vector
*			start - start index of long integral
*			stop - stop index of both integrals
*			lower - lower index range of short integral
*			upper - upper index range of short integral
*			method - see above...
*	----------------------------------------------------
*/
int PulseAnalysis::OptimizePSD (const std::vector< float > pulse, int start, int stop, int lower, int upper, float* paraL, float* paraS)
{
	if ((upper - lower) < 0) {return -1;}
	if ((upper - lower) < sizeof(paraS)/sizeof(float)) {return -1;}

	int j = 0;
	for (int i = 0; i < pulse.size(); i++) {
		if(pulse[i] > j) { j = (int)pulse[i]; l= i;}
	}

	// Constant fraction discrimination (CFD) at 50% of amplitude

	j = l;
	for(int i = j - 50; i < j; i++) {
		if(pulse[i] < (pulse[l])*0.5) {k = i;}
	}

	x = ((0.5*(pulse[l]) - pulse[k - 1])/((pulse[k + 1]
			- pulse[k - 1])/3))+ ((float)k - 1);


	if((k - start) > 0 && (k + start) < pulse.size()) {
			// Initialization
			integralL = 0;
			for ( j = 0; j < upper - lower; j++){
				 paraS[j] = 0;
			}

			// Begin integration using trapezoidal rule
			for (int i = (k - start); i < (k + stop); i++) {

				integralL += 0.5*(pulse[i-1] + pulse[i]);

				if (i >= (k + lower)) {
					for ( j = 0; j < upper - lower; j++)
					{
						if (i >= (k + lower + j)) {
							 paraS[j] += 0.5*(pulse[i-1] + pulse[i]);
						}
					}
				}

			}
	}

	if (integralL != 0) {
		*paraL = integralL;
		*paraS = integralS;
		return 0;
	} else { return -1;}

}


/** ----------------------------------------------------
*	Determine PSD by charge integration method
*		- Calculates discrimination parameters by
*		  integration of user defined regions of an
*		  input pulse.
*
*	Methods:
*			1 - trapezoidal rule
*			2 - composite Simpson's rule
*			3 - rectangular method
*			4 - Harwell method
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			start - start index of long integral
*			stop - stop index of both integrals
*			offset - start index of short integral
*			method - see above...
*	----------------------------------------------------
*/
int PulseAnalysis::PSD_Integration (const std::vector< float > pulse, int start, int stop, int offset, int method, float* paraL, float* paraS)
{
	int j = 0;
	for (size_t i = 0; i < pulse.size(); i++) {
		if(pulse[i] > j) { j = (int)pulse[i]; l= i;}
	}

	// Constant fraction discrimination (CFD) at 50% of amplitude

	j = l;
	for(int i = j - 50; i < j; i++) {
		if(pulse[i] < (pulse[l])*0.5) {k = i;}
	}

	x = ((0.5*(pulse[l]) - pulse[k - 1])/((pulse[k + 1]
			- pulse[k - 1])/3))+ ((float)k - 1);

	if((k - start) > 0 && (k + start) < pulse.size()) {
			// Initialization
			integralS = 0; integralL = 0;

			// Begin integration using trapezoidal rule
			if (method == 1) {
				for (int i = (k - start); i < (k + stop); i++) {

					integralL += 0.5*(pulse[i-1] + pulse[i]);

					if (i > (k + offset)) { integralS += 0.5*(pulse[i-1] + pulse[i]);}
				}
			}

			// Begin integration using composite Simpson's rule
			if (method == 2) {

				sumL = 0;
				sumS = 0;
				integralL = pulse[k - start];
				integralL = pulse[k + offset];

				for (int i = (k - start) + 1; i < (k + stop) - 2; i+=2) {

					sumL += pulse[i];

					if (i > (k + offset)) { sumS += pulse[i];}
				}

				integralL += 4*sumL;
				integralS += 4*sumS;

				sumL = 0;
				sumS = 0;
				for (int i = (k - start + 2); i < (k + stop) - 3; i+=2) {

					sumL += pulse[i];

					if (i > (k + offset)) { sumS += pulse[i];}
				}

				integralL += 2*sumL;
				integralS += 2*sumS;

				integralL += pulse[k + stop];
				integralS += pulse[k + stop];

				integralL = integralL/3;
				integralS = integralS/3;
			}

			// Begin integration using rectangular method
			if (method == 3) {
				for (int i = (k - start); i < (k + stop); i++) {

					integralL += pulse[i];

					if (i > (k + offset)) { integralS += pulse[i];}
					//if (i > (k - start) && i <= (k + offset)) { integralS += pulse[i];}
				}

				// Adjust for error due to array indexing (i.e. rounding)
				if (x < k + offset)
				{
					//integralS += ((pulse[k + offset] - pulse[k + offset-1])*x*0.5) + ((x - (float)(k + offset))*pulse[k + offset-1]);
				}

			}

	} else {return -1;}

	if (integralS != 0) {
		*paraL = integralL;
		*paraS = integralS;
		return 0;
	} else { return -1;}
}

/** ----------------------------------------------------
*	Restores the baseline of a input pulse
*
*	Methods:
*			1 - SNIP fitting routine
*			2 - baseline averaging
*
*	Inputs:
*			pulse - input array
*			iterations - number of iterations
*			method - see above...
*
*	----------------------------------------------------
*/
int PulseAnalysis::Baseline_restore (std::vector< float > & pulse, float* baseline, int iterations, int method)
{
	if (method == 1) {
		for (size_t i=0; i< pulse.size(); i++) {
			baseline[i] = pulse[i] + 0.3;
		}

		// Michael Febbraro 3/2/2013
		cout << "Method Disabled!" << endl; return -1;

		//TSpectrum::Background(baseline, length, iterations,TSpectrum::kBackDecreasingWindow,TSpectrum::kBackOrder2,
		//		TSpectrum::kTRUE, TSpectrum::kBackSmoothing3,TSpectrum::kFALSE);

		for (size_t i=0; i< pulse.size(); i++) {
			pulse[i] = pulse[i] - baseline[i];
		}
	}

	if (method == 2) {
		x = 0;
		for (int i = 0; i < iterations; i++) {x += pulse[i];}
		x = x / (float)iterations;

		for (size_t i = 0; i < pulse.size(); i++) {
			pulse[i] -= x;
		}
	}

	if (method == 3) {
		int j = 0;
		for (size_t i = 0; i < pulse.size(); i++) {
		for (auto value : pulse) {
			if(value > j) { j = (int)value; k = i;}
		}

		x = 0; y = 0;
		for (int i = 0; i < iterations; i++) {x += pulse[i];}
		x = x/(float)iterations;

		for (auto & value : pulse) {
			value -= x;
		}
	}

	return 0;
}

/** ----------------------------------------------------
*	Determine time characteristics of the a pulse
*		- Fits timing regions using linear regression
*		  over a user defined range.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			range - range of linear regression for timing
*	----------------------------------------------------
*/
int PulseAnalysis::Parameters (float* pulse, int length, int range, float* CFD, float* amplitude, float* risetime, float* falltime, float* width)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	int j = 0;
	for (int i = 0; i < length; i++) {
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}

	*amplitude = pulse[k];

	// Constant fraction discrimination (CFD) at 50% of amplitude
	memset(linear, 0, sizeof(linear));


	j = k;
	for(int i = j - 50; i < j; i++) {
		if(pulse[i] < (*amplitude)*0.5) {k = i;}
	}

	/*
	for(i = (point_t - (int)(range/2)); i <= (point_t - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;

	*CFD = ((0.5*(*amplitude) - pulse[(point_t - (int)(range/2))])/((linear[0]
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1])))+ (point_t - (int)(range/2));

	*/

	*CFD = ((0.5*(*amplitude) - pulse[(k - (int)(range/2))])/((pulse[k + (int)(range/2)] - pulse[k - (int)(range/2)])/(float)range))+ ((float)k - ((float)range/2));


	// Risetime by 10-90% Method
	j = k;
	for(int i = j - 50; i < j; i++)
	{
		if(pulse[i] < (*amplitude)*0.9) {k = i;}
	}

	memset(linear, 0, sizeof(linear));
	for(int i = (k - (int)(range/2)); i <= (k - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;

	*risetime = ((0.9*(*amplitude) - pulse[(k - (int)(range/2))])/((linear[0]
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1]))) + (k - (int)(range/2));

	j = k;
	for(int i = j - 50; i < j; i++) {
		if(pulse[i] < (*amplitude)*0.1) {k = i;}
	}

	memset(linear, 0, sizeof(linear));
	for(int i = (k - (int)(range/2)); i <= (k - (int)(range/2) + range); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/range;
	linear[1] = linear[1]/range;
	linear[2] = linear[2]/range;
	linear[3] = linear[3]/range;

	*risetime = *risetime - ((0.1*(*amplitude) - pulse[(k - (int)(range/2))])/((linear[0]
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1]))) + (k - (int)(range/2));


		// NOTE: Still need to do fall time and width...
		*falltime = 0;
		*width = 0;

		return 0;
}

/** ----------------------------------------------------
*	Determine time characteristics of the a pulse
*		- Fits timing regions using linear regression
*		  over a user defined range.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			range - range of linear regression for timing
*	----------------------------------------------------
*/
int PulseAnalysis::Parameters2 (float* pulse, int length, int range, float* CFD, float* amplitude)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	int j = 0;
	m=0;
	for (int i = 0; i < length; i++) {
		if(pulse[i] > j) { j = (int)pulse[i]; k = i; m = i;}
	}

	*amplitude = pulse[k];

	// Constant fraction discrimination (CFD) at 50% of amplitude
	float frac = 0.5;
	j = k;
	for(int i = j - 10; i < j; i++)
	  {
	    if(pulse[i] > (*amplitude)*0.5)
	      {
		*CFD = (pulse[m]*0.5 - pulse[i] + (pulse[i]-pulse[i-1])*(float)i ) / (pulse[i]-pulse[i-1]);
		break;
	      }
	  }

	// CFD by 10-90% Method
	//float x10, y10, x90, y90;
	/*
	j = k;
	for(i = j - 10; i < j; i++)
	  {
	    if(pulse[i] < (*amplitude)*0.1) {k = i;}
	  }

	x10 = ((pulse[m]*0.1 - pulse[k]) / (pulse[k+1] - pulse[k])) + (float)k;
	y10 = pulse[m]*0.1;

	j = m;
	for(i = j - 10; i < j; i++)
	  {
	    if(pulse[i] < (*amplitude)) {k = i;}
	  }

	x90 = ((pulse[m]*0.9 - pulse[k]) / (pulse[k+1] - pulse[k])) + (float)k;
	y90 = pulse[m]*0.9;

	*CFD = ( pulse[m]*frac - y10 + x10*((y90-y10)/(x90-x10)) ) / ((y90-y10)/(x90-x10));
	*/

	return 0;
}

/** ----------------------------------------------------
*	Determine location and number of peaks in a waveform
*		- Determines the number of peaks by applying a
*		  above-threshold condition at a user defined
*		  number of std. dev. of the baseline all using
*		  the first derv. of the input waveform.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			sigma - threshold in number of std. deviations
*			range - number of points to deterimine std. dev.
*	Outputs:
*			numPeaks - number of peaks found
*			locPeaks - array of array indexes of peak locations
*
*	----------------------------------------------------
*/
int PulseAnalysis::PeakFinder (float* pulse, int length, int sigma, int range, int method, int* numPeaks, int* value)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	//locPeaks = new int*[10];
	for(int i = 0; i < 10; ++i) {
		//locPeaks[i] = new int[3];
	}

// 	// Calculate first derivative and average of first derv.
	x = 0;
	for (int i = 0; i < length - 1; i++) {
		pulse_temp[i] = pulse[i + 1] - pulse[i];
		x += pulse_temp[i];
	}


	// Determine standard deviation of baseline
	x = x/range;
	z = 0;
	for (int i = 0; i < range; i++) { z += pow((pulse_temp[i] - x), 2); }
	z = sqrt(z/range);

	*numPeaks = 0;

	if(method == 1) { w = sigma*z; }

	if(method == 2) { w = sigma; }

	int j = 0;
	k = 0; x = 0; y = 0; m = 0; n = 0;
	positive = 0; negative = 0;
	for (int i = 0; i < length - 5; i++) {
		if( pulse_temp[i] <= w && pulse_temp[i] >= -1*w) {
			if (positive == 1 && negative == 1) {
				positive = 0;
				negative = 0;
				//locPeaks[*numPeaks - 1][1] = j;
				//locPeaks[*numPeaks - 1][2] = k;
				//locPeaks[*numPeaks - 1][0] = n;
				*value = k;

				x = 0; y = 0;
			}
		}
		else {
			if (positive == 0 && negative == 0 && i > n + 40) { m++; n = i;}

			// Positive inflection point
			if(pulse_temp[i] > 0) {
				positive = 1;
				if (pulse_temp[i] > x) {x = pulse_temp[i]; j = i;}
			}

			// Negative inflection point
			if(pulse_temp[i] < 0)
			{
				negative = 1;
				if (pulse_temp[i] < y) {y = pulse_temp[i]; k = i;}
			}

		}

	}

	*numPeaks = m;

	return 0;
}

/** ----------------------------------------------------
*	Determine the nth derivative of a waveform
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			order - order of the dervative (i.e. 1 = first, 2 = second, ...)
*	----------------------------------------------------
*/
int PulseAnalysis::Derivative(float* pulse, int length, int order)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	for(int j = 1; j <= order; j++) {
		for (int i = 0; i < length - j; i++) {
			pulse[i] = pulse[i + 1] - pulse[i];
		}
	}
	return 0;
}

/** ----------------------------------------------------
*	Determine the nth derivative of a waveform
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			order - order of the dervative (i.e. 1 = first, 2 = second, ...)
*	----------------------------------------------------
*/
int PulseAnalysis::Integral(float* pulse, int length)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	for (int i = 0; i < length - 1; i++)
	{
			pulse[i + 1] = pulse[i + 1] + pulse[i];
	}
	return 0;
}

/** ----------------------------------------------------
*	Output a time pickoff from user defined condition
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*
*	----------------------------------------------------
*/
int PulseAnalysis::Time_Pickoff (float* pulse, int length, int range, int low, int high , int method, float* CFD)
{
	float max, time;
	int index;

	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	if (method == 1)
	  {
	    max = 0;

	    for (int i = low; i <= high; i++)
	      {
		if (pulse[i] > max)
		  {
		    max = pulse[i];
		    index = i;
		  }
	      }

	// Constant fraction discrimination (CFD) at 50% of amplitude

	    for(int i = low; i < index; i++)
	      {
		if(pulse[i] < (max*0.5)) {k = i;}
	      }

	    time = (((max*0.5) - pulse[k])/(pulse[k+1] - pulse[k])) + (float)k;

	    if (time >= low && time <= high){*CFD = time;}
	    else {time = -1;}

	  }

	if (method == 2)
	  {

	    max = 0;

	// Find first value over 'range'
	    for (int i = low; i <= high; i++)
	      {
		if (pulse[i] > range)
		  {
		    k = i; break;
		  }
	      }

	    time = (((((float)(range)) - pulse[k-1]))/(pulse[k] - pulse[k-1])) + (float)k;

	    if (time >= low && time <= high){*CFD = time;}
	    else {time = -1;}


	  }

	return 0;
}

int PulseAnalysis::PSD_Zerocross (float* pulse, int length, int integration, int differentiation, float* PSD)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	int j = 0;
	for (int i = 0; i < length; i++) {
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}

	// CFD timing pickoff at 50% of amplitude
	memset(linear, 0, sizeof(linear));

	for(int i = (k - (int)(3/2)); i <= (k - (int)(3/2) + 3); i++) {
		linear[0] += i*pulse[i];
		linear[1] += i;
		linear[2] += pulse[i];
		linear[3] += i*i;
	}
	linear[0] = linear[0]/3;
	linear[1] = linear[1]/3;
	linear[2] = linear[2]/3;
	linear[3] = linear[3]/3;

	*PSD = ((0.5*(pulse[k]) - pulse[(k - (int)(3/2))])/((linear[0]
		- linear[1]*linear[2])/(linear[3] - linear[1]*linear[1])))+ (k - (int)(3/2));

	// Shaping amplifier
	memset(pulse_temp, 0, sizeof(pulse_temp));
	for(int i = (int)integration; i < (length - integration); i++) {
		for(j = 0; j < integration; j++) {
			if((j + i) >= 0 && (j + i) < length) {pulse_temp[i] = pulse_temp[i] + pulse[j + i];}
		}
	}

	//memset(pulse, 0, sizeof(pulse));
	//for (i = 0; i < (length - 1); i++) {pulse[i] = pulse_temp[i + 1] - pulse_temp[i];}


	//PSD = TMath::LocMin(length,pulse) - TMath::LocMax(length,pulse);


	memset(pulse_temp, 0, sizeof(pulse_temp));
	for (int i = 0; i < (length - 1); i++) {pulse_temp[i] = pulse[i + 1] - pulse[i];}

	//for (int i=0; i< length; i++) {pulse[i] = pulse_temp[i];}
	return 0;
}

/** ----------------------------------------------------
*	Determine the integral below half the pulse amplitude.  Used for
*   reconstruction of partial pulses.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			integral - integral below half the pulse amplitude
*	----------------------------------------------------
*/
int PulseAnalysis::Half_Integral(float* pulse, int length, float* integral)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	*integral = 0;
	int j = 0;
	for (int i = 0; i < length; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; k = i;}
	}

	for (int i = 0; i < length ; i++)
	{
		if(pulse[i] < (0.5*pulse[k]))
		{
			*integral = pulse[i + 1] + pulse[i];
		}
		else
		{
			*integral = pulse[i + 1] +pulse[k];
		}
	}


	return 0;
}

/** ----------------------------------------------------
*	Smoothing method for waveforms
*
*	Methods:
*			1 - Moving Average (option = range of moving average)
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			iterations - number of iterations to apply smoothing method
*			method - see above
*			option - see above
*	----------------------------------------------------
*/
int PulseAnalysis::Smooth(float* pulse, int length, int iterations, int method, float option)
{
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	for(k = 0; k < iterations; k++)
	{
	// Moving average
	if (method == 1)
	{
		for (int i = 2; i < length - 3; i++)
		{
			MovingAverage[i] = 0;
			for (int j = i - 2; j <= i + 2; j++) {MovingAverage[i] += pulse[j];}
		}
		for (int i =1; i < length - 3; i++)
		{
			pulse[i] = (MovingAverage[i]/5) ;
		}
	}

	}
	return 0;
}

/** ----------------------------------------------------
*	HPGe
*	- This method processes signals from HPGe detectors.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			amplitude - amplitude of the pulse
*	----------------------------------------------------
*/
int PulseAnalysis::HPGe(float* pulse, int length, float* amplitude){
	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	return 0;
}

/** ----------------------------------------------------
*	PSD routine for pulses with after pulsing
*		- Calculates discrimination parameters by
*		  integration of user defined regions of an
*		  input pulse. Note this method performs short
*         integral on front egde of pulse.
*
*	Inputs:
*			pulse - input array
*			length - the length of pulse
*			start - start index of long integral
*			stop - stop index of both integrals
*			offset - start index of short integral
*			apoffset - index for begining after pulsing threshold
*           apthres - threshold for after pulsing trigger
*	----------------------------------------------------
*/
int PulseAnalysis::PSD_Integration_Afterpulsing (float* pulse, int length, int start, int stop, int offset, int apoffset, float apthres, int* pposition, float* paraL, float* paraS, int* ap)
{

	if (length < sizeof(pulse)/sizeof(float)) {return -1;}

	int j = 0;
	for (int i = 0; i < length; i++)
	{
		if(pulse[i] > j) { j = (int)pulse[i]; l= i;}
	}

	// Constant fraction discrimination (CFD) at 50% of amplitude

	j = l;
	for(int i = j - 50; i < j; i++)
	{
		if(pulse[i] < (pulse[l])*0.5) {k = i;}
	}

	//*pposition = l;

	x = ((0.5*(pulse[l]) - pulse[k - 1])/((pulse[k + 1]
			- pulse[k - 1])/3))+ ((float)k - 1);


	if((k - start) > 0 && (k + start) < length) {
			// Initialization
			integralS = 0; integralL = 0; *ap = 0; j=0;

			// Begin integration using rectangular method
            for (int i = (k - start); i < (k + stop); i++) {

                integralL += pulse[i];

                if (i >= (k + offset)) { integralS += pulse[i];}
                if (i > (k + apoffset) && pulse[i] > j) {
                    j = pulse[i]; *ap = j;
                    *pposition = i;
                }
            }

	} else {return -1;}

	if (integralS != 0) {
		*paraL = integralL;
		*paraS = integralS;
		return 0;
	} else { return -1;}

}
