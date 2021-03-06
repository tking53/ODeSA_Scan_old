/************************************************************************
 *
 *  Filename: Scan.cpp
 *s
 *  Description:
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     David Walter
 *     Karl Smith
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 5/11/2021
 *
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include "PulseAnalysis.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"
#include "TF1.h"
#include <TGraph.h>

#include <cxxopts.hpp>

using namespace std;

typedef struct
{
  Float_t l;               // Long integral
  Float_t s;               // Short integral
  Float_t amp;             // Amplitude
  Float_t cfd;             // Trigger time
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Detector trigger

} Can;

enum class CfdMode {Lin, Fit};

float CfdFit(const std::vector< float > & pulse, int pposition);
float CfdLin(const std::vector< float > & pulse, int pposition);
double Rf_KS(const std::vector< float > & pulse, float thresh=-1);

int Scan (std::string prefix, std::string filename, std::string treeDesc="", CfdMode cfdMode=CfdMode::Lin, int num=0,
	std::vector< short > traceChannels = {}){
	const int numDets = 16;
	static std::array< Can, numDets > dets;

	std::array< ifstream, numDets > fps;

	long TEvt = 0;

	uint32_t buffer32;
	uint16_t buffer16;

	Double_t rf = -1;
	Double_t ebit = -1;

  TSpectrum *s = new TSpectrum();


  /** ----------------------------------------------------
   *	Calibrations and threshold
   *   ----------------------------------------------------
   */

  std::array< float, numDets > cal =
  {   0.0359,
      0.0327,
      0.0362,
      0.0269,
      0.0327,
      0.0345,
      0.0394,
      0.0376,
      0.031,
      0.0236,
      1.0, 1., 1., 1., 1., 1.

  }; // calibration (keVee / (bit * sample)) from manual check on calibration

  std::array< float, numDets > caloffset =
  {	21.124,
	23.803,
	18.688,
	18.465,
	11.81,
	23.183,
	22.302,
	18.863,
	13.401,
	30.319,
	0,
	0,
	0,
	0,
	0,
	0
  };

  std::array< float, numDets > threshold =
  {	15940,
	15900,
	15840,
	15930,
	15920,
	15970,
	15930,
	15970,
	15950,
	15850,
	15750,
	15850,
	158250,
	15700,
	15700,
	15700
  };

  /** ----------------------------------------------------
	*	Get functions
	*   ----------------------------------------------------
	*/

  PulseAnalysis *Analysis = new PulseAnalysis();


  TFile *ff = new TFile(filename.c_str(), "RECREATE");

  TTree *tt = new TTree("T", treeDesc.c_str());

	for (int i=0;i < dets.size(); i++) {
		tt->Branch(("d" + std::to_string(i)).c_str(), &dets.at(i),"l:s:amp:cfd:psd:trg");
	}

	std::array< TH1F*, numDets > traceHists{nullptr};
	for (auto ch : traceChannels) {
		std::cout << "Storing traces for channel: " << ch << std::endl;
		std::string name = "trace" + std::to_string(ch);
		TH1F *trace = new TH1F(name.c_str(),
		                       ("Trace for channel " + std::to_string(ch)).c_str(), 200, 0, 199);
		traceHists[ch] = trace;
		tt->Branch(name.c_str(), "TH1F", &traceHists[ch]);
	}

  //tt->Branch("runtime",&runtime,"Runtime (ms)");     // Runtime in ms
  tt->Branch("rf",&rf, "rf/D"); // RF timing signal
  tt->Branch("ebit",&ebit, "ebit/D"); // RF timing signal

  int numSteps = 20;
  float psd_opt[20];
  //tt->Branch("psd_opt", &psd_opt, "psd_opt[20]/F");

	// Open files
	bool valid = false;
	size_t file_length = 0;
	for (int i = 0; i < numDets; i++) {
		std::string openfile = prefix + "_wave" + to_string(i) + ".dat";
		cout << " Opening file: " << openfile;
		fps[i].open(openfile, std::ifstream::in | std::ifstream::binary);

		if (fps[i].is_open()) {
			cout << " - Open!" << endl;
			valid = true;
		}
		else {cout << " - Failed!" << endl;}
		if (! file_length && fps[i].is_open()) {
			fps[i].seekg(0, std::ios::end);
			file_length = fps[i].tellg();
			fps[i].seekg(0, std::ios::beg);
		}
	}
	if (! valid) {
		cerr << "No valid files found!" << endl;
		exit(-1);
	}

	float runtime = 0;
	float prevtime = 0;
	float trgtime = 0;
	float difftime = 0;

	/** ----------------------------------------------------
	 *	Process liquid scintillator det events
	 *   ----------------------------------------------------
	 */

	bool data_valid = true;
	while (data_valid) {
		// Stop after nth events...
		if (num && TEvt > num) {break;}

		rf = -1;
		ebit = -1;

		// Loop over each channel
		for (int detNum=0; detNum < numDets; detNum++) {

		  // Binary parsing
			auto & fp = fps[detNum];

			if ( ! fp.is_open()) {continue;}
			if (!fp.read((char*)&buffer32, 4)) {
				data_valid = false;
				break;
			}
			int Tracelength = (buffer32 - 16)/2;
			fp.read((char*)&buffer32, 4);
			fp.read((char*)&buffer32, 4);
			fp.read((char*)&buffer32, 4);

			// Trigger time in 2 ADC clock cycles ( 8 ns )
			if (detNum == 0) trgtime = buffer32;

			// Reset variables
			float CFD = -1;
			float amplitude = -1;
			float paraL = 0;
			float paraS = 0;
			bool trg = false;
			int pposition = -1;

			std::vector< float > pulse;
			pulse.reserve(Tracelength);

			// Get traces
			for (int i = 0; i < Tracelength; i++) {
				if (!fp.read((char*)&buffer16, 2)) {
					data_valid = false;
					break;
				}

				// Determine the trigger channel by inspecting  
				if (buffer16 > threshold[detNum]) {trg = true;}

				// Flip detector signals to positive
				if (detNum < 12) {
					pulse.push_back(16383 - (float) buffer16);
				}
				else {pulse.push_back(buffer16);}
			}

			// If a histogram is defined for the trace we set the values.
			if (traceHists[detNum]) {
				TH1* trace = traceHists[detNum];
				for (int i = 0; i < Tracelength; i++) {
					trace->SetBinContent(i, pulse[i]);
				}
			}


			/** Liquid can processing **/
			if (Tracelength > 1) {
				if (detNum < 13) {
					// Process trace
					// Get a Continuous Moving Average (CMA) of the pulse.
					std::vector< float > CMAtrace = Analysis->CMA_Filter(pulse, 10, pulse[0], 3.5 );
					// Subtract the CMA from the pulse
					std::transform(pulse.begin(), pulse.end(), CMAtrace.begin(), pulse.begin(), std::minus<int>());
				}

				// Find the maximum element in the subtracted pulse.
				auto maxElement = std::max_element(pulse.begin(), pulse.end());
				amplitude = *maxElement;
				pposition = std::distance(pulse.begin(), maxElement);

				if (trg) {
					if (detNum < 13) {
						switch (cfdMode) {
							case CfdMode::Fit:
								CFD = CfdFit(pulse, pposition);
								break;
							case CfdMode::Lin:
								CFD = CfdLin(pulse, pposition);
								break;
						}
					}
					else {
						CFD = (float)pposition;
					}
				}
				else {CFD = -1;}

				// PSD integration
				float offset = 12.0; // original 12
				if (pposition - 10 > 0 && pposition + 100 < Tracelength) {
					for (int i = (pposition - 10); i < (pposition + 100); i++) {
						paraL += pulse[i];
						if (i > pposition + offset) { paraS += pulse[i];}
					}
				}
			}

			auto & det = dets.at(detNum);
			det.amp = amplitude;
			det.l = paraL * cal[detNum];
			det.s = paraS;
			det.cfd = CFD;
			det.trg = trg;
			det.psd = paraS / paraL;
			if (det.trg) {det.trg = 1;}
			else {det.trg = 0;}

			switch(detNum) {
				case 0 :
					// Runtime clock
					difftime = trgtime - prevtime;
					if(difftime < 0) {
						difftime = difftime + 2147483647;
					runtime += ((8*difftime)/1.0E6);
					prevtime = trgtime;
					}
					else {
						runtime += ((8*difftime)/1.0E6);
						prevtime = trgtime;
					}
					break;

				case 14:
					rf = Rf_KS(pulse);
					break;

				case 15:
					ebit = pulse[0];
					break;

				default:
					break;
			}

		}
		tt->Fill();
      TEvt++;
      if (TEvt % 1000 == 0) {
			int percentComplete = fps[1].tellg() / (float) file_length * 100;
            // To fix - line 181 make file_length size of biggest file not last file
			cout << "\rEvent: " << TEvt << " " << percentComplete << "%" << flush;
		}

	}

	for (auto & fp : fps) {
		if(fp.is_open()) fp.close();
	}
	ff->cd();
	tt->Write();
	ff->Close();

	cout << "\nFinished! " << TEvt << " events" << endl;

	return 0;
}

double Rf_BT(const std::vector<float> & pulse) {
	for(size_t i=0; i<pulse.size(); i++){
		if(pulse[i]-pulse[i-1] < 200) continue;
		// ibt-1 to ibt+2 linear fit
		double slope = (pulse[i+2] - pulse[i-1]) / 3.;
		double intercept = pulse[i+2] - (slope * ((float)i+2));
		double thresh = ((pulse[i+2]-pulse[i-1])/2.) + pulse[i-1];
		printf("%f %f %f\n", slope, intercept, thresh);
		return (thresh - intercept) / slope;
	}
	return 0;
}

/**Computes the x-position where an RF sine signal crosses a threshold.
 * If not provided the threshold is computed as the midway point bewteen the
 * max and min. The provided signal is scanned until the values are below the
 * threshold and then the search for the crossing begins. This ensures that the
 * crossing is on a rising value. The crossing point is then linearly
 * interpolated between the two bracketing values.
 *
 * \param[in] pulse The signal to scan.
 * \param[in] thresh An optional threshold to determine the crossing. If less
 * 	than 0 the value is computed form the pulse.
 * \return The x-position of the threshold crossing.
 */
double Rf_KS(const std::vector< float > & pulse, float thresh /* = -1 */) {
	// Compute a threshold if not provided.
	if (thresh < 0) {
		float max = *std::max_element(pulse.begin(), pulse.end());
		float min = *std::min_element(pulse.begin(), pulse.end());
		thresh = (max - min) / 2. + min;
	}

	// Find a point below the threshold to start scan for a threshold crossing.
	size_t start;
	for (size_t i=0; i<pulse.size(); ++i) {
		if (pulse[i] < thresh) {
			start = i;
			break;
		}
	}
	// Find the threshold crossing
	for (size_t i=start; i<pulse.size(); ++i) {
		if (pulse[i] < thresh) continue;
		double slope = (pulse[i] - pulse[i-1]) / 2.;
		double intercept = pulse[i] - (slope * i);
		return (thresh - intercept) / slope;
	}

	return -1;

}

float CfdFit(const std::vector< float > & pulse, int pposition) {
	TGraph g;
   int start = pposition - 5;
   int stop = pposition + 2;
	for (size_t i=start; i < stop; ++i) g.SetPoint(i - start, i, pulse[i]);

	static TF1 *f1 = new TF1("f1", "gaus", 0, pulse.size());
	f1->SetParameters(1.0, (double)pposition, 0.1);

	g.Fit("f1","Q");

	float mu = (float)f1->GetParameter(1);
	float sigma = (float)f1->GetParameter(2);

	return mu - sqrtf(1.38629*sigma*sigma);
}

float CfdLin(const std::vector< float > & pulse, int pposition) {
	// Timing signal using a linear interpolation
	float xy_ = 0;
	float x_ = 0;
	float y_ = 0;
	float x_2 = 0;

	for(int i = pposition - 3; i <= pposition + 3; i++) {
		x_ += i;
		x_2 += i * i;
		y_ += pulse[i];
		xy_ += pulse[i] * i;
	}

	float m_ = ((7.0*xy_) - (x_*y_))/((7.0*x_2) - (x_*x_));
	float b_ = ((y_*x_2) - (x_*xy_))/((7.0*x_2) - (x_*x_));

	return (0.5 * pulse[pposition] - b_) / m_;
}


int main(int argc, char ** argv) {
	std::string filename, treeDesc, prefix;
	std::vector< short > traceChannels;

	//Parser command arguments
	cxxopts::Options parser(argv[0], std::string("E20003 Scan"));
	parser.add_options("")
		("n,num", "Number of events", cxxopts::value< int >()->default_value("0"), "NUM")
		("desc", "Tree description", cxxopts::value(treeDesc)->default_value(""), "DESC")
		("cfd", "CFD Mode (lin, fit)", cxxopts::value< std::string >()->default_value("lin"), "MODE")
		("trace", "Store traces", cxxopts::value(traceChannels), "CH")
		("i,interactive", "Run interactively")
		("h,help", "Help");

	// Positional arguments
	parser.parse_positional({"prefix", "output"});
	parser.positional_help("<prefix> <output>");
	parser.add_options("positional")
		("prefix", "Prefix", cxxopts::value(prefix)->default_value(""), "PREFIX")
		("output", "Output", cxxopts::value(filename)->default_value(""), "OUTPUT");

	// parse arguments
	auto args = parser.parse(argc, argv);

	// Print help message if requested
	if (args["help"].as<bool>()) {
		std::cerr << parser.help({""});
		return 0;
	}


	CfdMode cfdMode;
	if (args["cfd"].as<std::string>() == "lin") {
		cfdMode = CfdMode::Lin;
	}
	else if (args["cfd"].as<std::string>() == "fit") {
		cfdMode = CfdMode::Fit;
	}
	else {
		std::cerr << "ERROR: Unknown CFD Mode: '" << args["cfd"].as<std::string>() << "'!" << std::endl;
		return 2;
	}

	/** ----------------------------------------------------
	 *	Program start...
	 *   ----------------------------------------------------
	 */

	cout << " ------------------------------------------------ " << endl;
	cout << " | Scan.cpp - binary version                     |" << endl;
	cout << " |   Experiment: 20Ne(d,n)                       |" << endl;
	cout << " |   Date: May 2021                              |" << endl;
	cout << " |   Calibration used: Becca 137Cs Pre BeamTime  |" << endl;
	cout << " |   ORNL Nuclear Astrophysics                   |" << endl;
	cout << " ------------------------------------------------ " << endl;

	if (args["interactive"].as<bool>()) {
		cout << "Root file name to be created: ";
		cin >> filename;

		cout << "Root tree description: ";
		cin >> treeDesc;

		cout << "Run binary file prefix: ";
		cin >> prefix;
	}

	if (filename == "") {
		std::cerr << "ERROR: An output file was not provided!" << std::endl;
		return 1;
	}

	return Scan(prefix,
	            filename,
	            treeDesc,
					cfdMode,
	            args["num"].as<int>(),
					traceChannels);
}
