#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <numeric>

using namespace std;
using namespace TMath;

class Hist2D
{
	private:
		size_t nxbins_, nybins_;
		double xmin_, xmax_;
		double ymin_, ymax_;
		vector<double> data_n;
		vector<double> data_w;
		vector<double> data_w2;
		vector<double> data_x;
		vector<double> data_y;
	public:
		Hist2D(size_t nxbins, double xmin, double xmax,
           size_t nybins, double ymin, double ymax)
           :
           nxbins_(nxbins), nybins_(nybins), xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
           data_n(nxbins*nybins, 0.), data_w(nxbins*nybins, 0.), data_w2(nxbins*nybins, 0.),
           data_x(nxbins*nybins, 0.), data_y(nxbins*nybins, 0.)
	  {
	  }

		int getXBin(double x) const {return (x-xmin_)*nxbins_/(xmax_ - xmin_);}
		int getYBin(double y) const {return (y-ymin_)*nybins_/(ymax_ - ymin_);}
		int getBin(double x, double y){return getXBin(x) + nxbins_ * getYBin(y);}
		void fill(double x, double y, double w = 1.)
		{
			int bin = getBin(x, y);
			if(bin < 0 || bin >= int(data_n.size())) return;

			data_n[bin] += 1;
			data_w[bin] += w;
			data_w2[bin] += w*w;
			data_x[bin] += x;
			data_y[bin] += y;
		}

		int getBinByBins(int xbin, int ybin) const {return xbin + nxbins_ * ybin;}

		double getBinEntries(int bin)
		{
			return data_n[bin];
		}
		double getBinEntries(int xbin, int ybin)
		{
			return getBinEntries(getBinByBins(xbin, ybin));
		}
		double getBinValue(int bin)
		{
			return data_w[bin];
		}
		double getBinValue(int xbin, int ybin)
		{
			return getBinValue(getBinByBins(xbin, ybin));
		}
		double getBinUnc(int bin)
		{
			if(data_n[bin] > 1)
			{
				return Sqrt(1./(data_n[bin] -1) *(data_w2[bin] - data_w[bin]*data_w[bin]/data_n[bin]));
			}
			else
			{
				return 0.;
			}
		}
		double getBinUnc(int xbin, int ybin)
		{
			return getBinUnc(getBinByBins(xbin, ybin));
		}
		double getBinMeanX(int bin)
		{
			if(data_n[bin] > 0)
			{
				return data_x[bin]/data_n[bin];
			}
			else
			{
				return (xmax_ - xmin_)/nxbins_ * (0.5 + (bin % nxbins_))+ xmin_; 
			}
		}
		double getBinMeanX(int xbin, int ybin)
		{
			return getBinMeanX(getBinByBins(xbin, ybin));
		}
		double getBinMeanY(int bin)
		{
			if(data_n[bin] > 0)
			{
				return data_y[bin]/data_n[bin];
			}
			else
			{
				return (ymax_ - ymin_)/nybins_ * (0.5 + (bin / nxbins_)) + ymin_; 
			}
		}
		double getBinMeanY(int xbin, int ybin)
		{
			return getBinMeanY(getBinByBins(xbin, ybin));
		}
		int getNumBinsX() const {return nxbins_;}
		int getNumBinsY() const {return nybins_;}
		int getNumBins() const {return nxbins_*nybins_;}
		TH2D* writeTH2D(string name)
		{
			TH2D* h_n = new TH2D((name+"_n").c_str(), (name+"_n").c_str(), getNumBinsX(), xmin_, xmax_, getNumBinsY(), ymin_, ymax_); 
			TH2D* h_w = new TH2D((name+"_w").c_str(), (name+"_w").c_str(), getNumBinsX(), xmin_, xmax_, getNumBinsY(), ymin_, ymax_); 
			TH2D* h_w2 = new TH2D((name+"_w2").c_str(), (name+"_w2").c_str(), getNumBinsX(), xmin_, xmax_, getNumBinsY(), ymin_, ymax_); 
			TH2D* h_x = new TH2D((name+"_x").c_str(), (name+"_x").c_str(), getNumBinsX(), xmin_, xmax_, getNumBinsY(), ymin_, ymax_); 
			TH2D* h_y = new TH2D((name+"_y").c_str(), (name+"_y").c_str(), getNumBinsX(), xmin_, xmax_, getNumBinsY(), ymin_, ymax_); 
			for(int bx = 0 ; bx < getNumBinsX() ; ++bx)
			{
				for(int by = 0 ; by < getNumBinsY() ; ++by)
				{
					h_w->SetBinContent(bx+1, by+1, data_w[getBinByBins(bx, by)]);
					h_w->SetBinError(bx+1, by + 1, getBinUnc(bx, by));
					h_w2->SetBinContent(bx+1, by+1, data_w2[getBinByBins(bx, by)]);
					h_n->SetBinContent(bx+1, by+1, data_n[getBinByBins(bx, by)]);
					h_n->SetBinError(bx+1, by+1, Sqrt(data_n[getBinByBins(bx, by)]));
					h_x->SetBinContent(bx+1, by+1, data_x[getBinByBins(bx, by)]);
					h_y->SetBinContent(bx+1, by+1, data_y[getBinByBins(bx, by)]);
				}
			}
			return h_w;
		}
		Hist2D(string filename, string name)
		{
			TDirectory* cdir = gDirectory;
			TFile* fin = TFile::Open(filename.c_str());
			TH2D* h_n = dynamic_cast<TH2D*>(fin->Get((name+"_n").c_str()));
			TH2D* h_w = dynamic_cast<TH2D*>(fin->Get((name+"_w").c_str()));
			TH2D* h_w2 = dynamic_cast<TH2D*>(fin->Get((name+"_w2").c_str()));
			TH2D* h_x = dynamic_cast<TH2D*>(fin->Get((name+"_x").c_str()));
			TH2D* h_y = dynamic_cast<TH2D*>(fin->Get((name+"_y").c_str()));
			nxbins_ = h_n->GetNbinsX();
			nybins_ = h_n->GetNbinsY();
			xmin_ = h_n->GetXaxis()->GetXmin();
			xmax_ = h_n->GetXaxis()->GetXmax();
			ymin_ = h_n->GetYaxis()->GetXmin();
			ymax_ = h_n->GetYaxis()->GetXmax();
			data_n.resize(nxbins_*nybins_);
			data_w.resize(nxbins_*nybins_);
			data_w2.resize(nxbins_*nybins_);
			data_x.resize(nxbins_*nybins_);
			data_y.resize(nxbins_*nybins_);
			for(int bx = 0 ; bx < getNumBinsX() ; ++bx)
			{
				for(int by = 0 ; by < getNumBinsY() ; ++by)
				{
					data_n[getBinByBins(bx, by)] = h_n->GetBinContent(bx+1, by+1);
					data_w[getBinByBins(bx, by)] = h_w->GetBinContent(bx+1, by+1);
					data_w2[getBinByBins(bx, by)] = h_w2->GetBinContent(bx+1, by+1);
					data_x[getBinByBins(bx, by)] = h_x->GetBinContent(bx+1, by+1);
					data_y[getBinByBins(bx, by)] = h_y->GetBinContent(bx+1, by+1);
				}
			}
			delete fin;
			cdir->cd();
		}

};

class Hist1D
{
	private:
		int nbins_;
		double min_, max_;
		vector<double> data_n;
		vector<double> data_w;
		vector<double> data_w2;
		vector<double> data_x;
	public:
		Hist1D(size_t nbins, double min, double max) : nbins_(nbins), min_(min), max_(max), data_n(nbins, 0.), data_w(nbins, 0.), data_w2(nbins, 0.), data_x(nbins, 0.)
	{
	}

		void fill(double x, double w = 1.)
		{
			int bin = (x-min_)*nbins_/(max_ - min_);
			if(bin < 0 || bin >= nbins_) return;

			data_n[bin] += 1;
			data_w[bin] += w;
			data_w2[bin] += w*w;
			data_x[bin] += x;
		}
		int getNumBins() const {return nbins_;}
		double getXMin() const {return min_;}
		double getXMax() const {return max_;}

		double getBinEntries(int bin)
		{
			return data_n[bin];
		}
		double getBinValue(int bin)
		{
			return data_w[bin];
		}
		double getBinUnc(int bin)
		{
			if(data_n[bin] > 1)
			{
				return Sqrt(1./(data_n[bin] -1) *(data_w2[bin] - data_w[bin]*data_w[bin]/data_n[bin]));
			}
			else
			{
				return 0.;
			}
		}
		double getBinMeanX(int bin)
		{
			if(data_n[bin] > 0)
			{
				return data_x[bin]/data_n[bin];
			}
			else
			{
				return (max_ - min_)/nbins_ * (0.5 + (bin % nbins_))+ min_; 
			}
		}
		
		TH1D* writeTH1D(string name)
		{
			TH1D* h_n = new TH1D((name+"_n").c_str(), (name+"_n").c_str(), getNumBins(), min_, max_); 
			TH1D* h_w = new TH1D((name+"_w").c_str(), (name+"_w").c_str(), getNumBins(), min_, max_); 
			TH1D* h_w2 = new TH1D((name+"_w2").c_str(), (name+"_w2").c_str(), getNumBins(), min_, max_); 
			TH1D* h_x = new TH1D((name+"_x").c_str(), (name+"_x").c_str(), getNumBins(), min_, max_); 
			for(int b = 0 ; b < getNumBins() ; ++b)
			{
				h_n->SetBinContent(b+1, data_n[b]);
				h_n->SetBinError(b+1, Sqrt(data_n[b]));
				h_w->SetBinContent(b+1, data_w[b]);
				h_w->SetBinError(b+1, getBinUnc(b));
				h_w2->SetBinContent(b+1, data_w2[b]);
				h_x->SetBinContent(b+1, data_x[b]);
			}
			return h_w;
		}

		Hist1D(string filename, string name)
		{
			TDirectory* cdir = gDirectory;
			TFile* fin = TFile::Open(filename.c_str());
			TH1D* h_n = dynamic_cast<TH1D*>(fin->Get((name+"_n").c_str()));
			TH1D* h_w = dynamic_cast<TH1D*>(fin->Get((name+"_w").c_str()));
			TH1D* h_w2 = dynamic_cast<TH1D*>(fin->Get((name+"_w2").c_str()));
			TH1D* h_x = dynamic_cast<TH1D*>(fin->Get((name+"_x").c_str()));
			nbins_ = h_n->GetNbinsX();
			min_ = h_n->GetXaxis()->GetXmin();
			max_ = h_n->GetXaxis()->GetXmax();
			data_n.resize(nbins_);
			data_w.resize(nbins_);
			data_w2.resize(nbins_);
			data_x.resize(nbins_);
			for(int bx = 0 ; bx < getNumBins() ; ++bx)
			{
					data_n[bx] = h_n->GetBinContent(bx+1);
					data_w[bx] = h_w->GetBinContent(bx+1);
					data_w2[bx] = h_w2->GetBinContent(bx+1);
					data_x[bx] = h_x->GetBinContent(bx+1);
			}
			delete fin;
			cdir->cd();
		}
		double integral() const {return accumulate(data_w.begin(), data_w.end(), 0.);}
		double getEntries() const {return accumulate(data_n.begin(), data_n.end(), 0.);}
		double scale(double val)
		{
			double res = 0.;
			for(double& v : data_w) {v*=val; res+=v;}
			for(double& v : data_w2) {v*=val*val;}
			return res;
		}

};
