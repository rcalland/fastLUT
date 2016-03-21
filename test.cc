#include <iostream>
#include "fastLUT3D.h"

#include "TCanvas.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TF1.h"

int main()
{
  const bool morton_mode = false;
  
  TRandom3 *rnd = new TRandom3(0);
  
  const int s = 128;

  TH3F *h_lut = new TH3F("h_lut","",s,0,100,s,0,100,s,0,100);
  fastLUT3D *lut = new fastLUT3D(s);
  lut->setBounds(0,100,0,100,0,100);
  
  int c = 0;

  TF1 *fb = new TF1("fb","gaus(0)",0,100);
  fb->SetParameters(1.0, 50, 25);
  
  for(uint32_t i = 0; i < s; ++i)
    for(uint32_t j = 0; j < s; ++j)
      for(uint32_t k = 0; k < s; ++k)
	{
	  //float rndm = (float)rnd->Rndm() * 100;
	  float rndm = fb->Eval(float(i * float(100./s)));
	  rndm += fb->Eval(float(j * float(100./s)));
	  rndm += fb->Eval(float(k * float(100./s)));

	  //std::cout << i << " " << j << " " << k << " " << rndm << std::endl;
	  
	  lut->setElement(i,j,k, rndm, morton_mode);
	  h_lut->SetBinContent(i+1, j+1, k+1, rndm);
	  //std::cout << lut->getElement(i, j, k, morton_mode) << std::endl;
	  c++;
	}
  
  TH3I *hist = new TH3I("morton","",s,0,s,s,0,s,s,0,s);
  
  lut->makeMortonPlot(*hist);
  TCanvas *c1 = new TCanvas();
  hist->Draw();
  c1->SaveAs("morton.C");

  // lambda magic!
  auto genRnd = [=]() { return (float) (5.0f + rnd->Rndm() * 90); };

  auto trials = 1000;
  auto N = 33000;
  TH1F *time_hist = new TH1F("th","",100,0,5);
  
  for (auto t = 0; t < trials; ++t)
    {
      TStopwatch clock;
      clock.Start();

      float intrp[3];
      for (auto i = 0; i < N; ++i)
	{
	  for (int r = 0; r < 3; ++r) intrp[r] = genRnd();
	  float val = lut->Interpolate( intrp[0], intrp[1], intrp[2], morton_mode );
	  //float val_lut = h_lut->Interpolate( intrp[0], intrp[1], intrp[2]);
	  
	  //std::cout << val << " " << val_lut << std::endl;
	}
      
      clock.Stop();

      time_hist->Fill(clock.RealTime() * 1000.0);
      std::cout << "Time for " << N << " interpolations was " << clock.RealTime() * 1000.0 << " ms" << std::endl;
    }

  TCanvas *c2 = new TCanvas();
  time_hist->Draw();
  c2->SaveAs("timetrial.C");

  delete lut;
  delete rnd;
}
