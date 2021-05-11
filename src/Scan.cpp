/************************************************************************
 *
 *  Filename: Scan.cpp
 *s
 *  Description:
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     David Walter
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 8/24/2018
 *
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum NewScan.cpp
 *      - if errors copiling in Mac OSX
 *        - remove -03 option
 *        - clang++ instead of g++
 *        - $(root-config --cflags --libs) instead of 'root-config --cflags --libs'
 *
 *
 *  If "error while loading shared libraries: libcore.so: ..." occurs, type
 *  "source `root-config --prefix`/bin/thisroot.sh"
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

int Scan (){
	static std::array< Can, 13 > dets;
	bool beamON, trg;

	const int numFiles = 13;
	std::array< ifstream, numFiles > fps;

	string fileheader;

	int k, pposition,
    Tracelength,
    eventlength;

	std::vector< float > pulse, CMAtrace;

  Float_t amplitude,
    risetime,
    falltime,
    width,
    CFD,
    tac,
    paraL,
    paraS,
    runtime;

	std::string filename;
  char prompt[10],
    runnum[250],
    interrputPrompt;

	std::string prefix, openfile;

  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

  uint32_t buffer32;
  uint16_t buffer16;

  TSpectrum *s = new TSpectrum();

  TF1 *f1 = new TF1("f1","gaus",0,300);

  /** ----------------------------------------------------
   *	Calibrations and threshold
   *   ----------------------------------------------------
   */

  std::array< float, 16 > cal =
  {  0.0359,
	  0.0327,
	  0.0362,
	  0.0269,
	  0.0327,
	  0.0345,
	  0.0394,
	  0.0376,
	  0.031,
	  0.0236,
	  1.0,
	  1.,
	  1.,
	  1.,
	  1.,
	  1.
  }; // calibration (keVee / (bit * sample)) from manual check on calibration

  std::array< float, 16 > caloffset =
  {  21.124,
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

  std::array< float, 16 > threshold =
  {  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700,
	  15700
  };

  /** ----------------------------------------------------
	*	Get functions
	*   ----------------------------------------------------
	*/

  PulseAnalysis *Analysis = new PulseAnalysis();


  /** ----------------------------------------------------
	*	Program start...
	*   ----------------------------------------------------
	*/

  cout << " ------------------------------------------------ " << endl;
  cout << " | Scan.cpp - binary version                     |" << endl;
  cout << " |   Experiment: 20Ne(d,n)                       |" << endl;
  cout << " |   Date: May 2021                              |" << endl;
  cout << " |   Calibration used: None                      |" << endl;
  cout << " |   ORNL Nuclear Astrophysics                   |" << endl;
  cout << " ------------------------------------------------ " << endl;

  cout << "Root file name to be created: ";
  cin >> filename;

  cout << "Root file header: ";
  cin >> fileheader;

  cout << "Run binary file prefix ('../run#'): ";
  cin >> prefix;

  TFile *ff = new TFile((filename + ".root").c_str(), "RECREATE");

  TTree *tt = new TTree("T", fileheader.c_str());

  for (int i=0;i < dets.size(); i++) {
	  tt->Branch(("d" + std::to_string(i)).c_str(), &dets.at(i),"l:s:amp:cfd:psd:trg");
  }
  //tt->Branch("runtime",&runtime,"Runtime (ms)");     // Runtime in ms

  TH1F *trace0 = new TH1F("trace0","Trace for channel 0",200,0,199);
  tt->Branch("trace0","TH1F", &trace0);
  TH1F *trace1 = new TH1F("trace1","Trace for channel 1",200,0,199);
  //tt->Branch("trace1","TH1F", &trace1);

  TH1F *traceCMA = new TH1F("traceCMA","Trace for CMA",200,0,199);
  //tt->Branch("traceCMA","TH1F", &traceCMA);

  TH1F *trace0C = new TH1F("trace0C","Corrected trace for channel 0",200,0,199);
  //tt->Branch("trace0C","TH1F", &trace0C);

  int numSteps = 20;
  float psd_opt[20];
  //tt->Branch("psd_opt", &psd_opt, "psd_opt[20]/F");

  // Open files
  for (int i = 0; i < numFiles; i++) {
	  openfile = prefix + "_wave" + to_string(i) + ".dat";
	  cout << " Opening file: " << openfile;
	  fps[i].open(openfile, std::ifstream::in | std::ifstream::binary);

	  if (fps[i].is_open()) {cout << " - Open!" << endl;}
	  else {{cout << " - Failed!" << endl;} }
  }

  runtime = 0;
  prevtime = 0;

  /** ----------------------------------------------------
	*	Process liquid scintillator det events
	*   ----------------------------------------------------
	*/

  while (true) {
	  beamON = 0;
	  for (int detNum=0; detNum < numFiles; detNum++) {
		  // Stop after nth events...
		  //if (TEvt > 1000000) {break;}

		  // Binary parsing
			auto & fp = fps[detNum];

			if (!fp.read((char*)&buffer32, 4)) {break;}
			Tracelength = (buffer32 - 16)/2;
			fp.read((char*)&buffer32, 4);
			fp.read((char*)&buffer32, 4);
			fp.read((char*)&buffer32, 4);

			// Trigger time in 2 ADC clock cycles ( 8 ns )
			if (detNum == 0) trgtime = buffer32;

			// Reset variables
			CFD = -1;
			amplitude = -1;
			paraL = 0;
			paraS = 0;
			trg = 0;
			tac = 0;
			pposition = -1;

			// Get traces
			for (int i = 0; i < Tracelength; i++) {
				if (!fp.read((char*)&buffer16, 2)) {break;}

				// Flip detector signals to positive
				if (detNum < 12) {
					pulse.push_back(16383 - (float) buffer16);
				}
				else {pulse.push_back(buffer16);}

				if (pulse.back() > (16383 - threshold[detNum])) {trg = 1;}

				// Added traces
				if (detNum==0) {trace0->SetBinContent(i, pulse.back());}
				else if (detNum==1) {trace1->SetBinContent(i, pulse.back());}
			}

			/** Liquid can processing **/
			if (Tracelength > 1) {
				// Process trace
				// Get a Continuous Moving Average (CMA) of the pulse.
				CMAtrace = Analysis->CMA_Filter(pulse, 10, pulse[0], 3.5 );
				// Subtract the CMA from the pulse
				std::transform(pulse.begin(), pulse.end(), CMAtrace.begin(), pulse.begin(), std::minus<int>());

				// Find the maximum element in the subtracted pulse.
				auto maxElement = std::max_element(pulse.begin(), pulse.end());
				amplitude = *maxElement;
				pposition = std::distance(pulse.begin(), maxElement);

				// Fill the bins.
				for (size_t i=0; i < pulse.size(); i++) {
					auto & value = pulse[i];
					if (detNum==0) {trace0->SetBinContent(i, value);}
					else if (detNum==1) {trace1->SetBinContent(i, value);}
					trace0C->SetBinContent(i, value);
				}
				for (size_t i=0; i < CMAtrace.size(); i++) {
					traceCMA->SetBinContent(i, CMAtrace[i]);
				}


				// CFD timing
				//f1->SetParameters(1.0, (double)pposition, 0.1);
				//trace0->GetXaxis()->SetRangeUser(pposition - 5, pposition + 1);
				//trace0->Fit("f1","RQ");
				//float mu = (float)f1->GetParameter(1);
				//float sigma = (float)f1->GetParameter(2);

				//CFD = mu - sqrtf(1.38629*sigma*sigma);
				//CFD = (float)pposition;

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
			det.l = paraL;
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

				case 13:
				case 14:
				case 15:

				default:
					break;
			}

		}
		tt->Fill();
      TEvt++;
      if (TEvt % 1000==0) {cout << "\rEvent counter: " << TEvt << flush;}

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

int main() {

  return Scan();
}
