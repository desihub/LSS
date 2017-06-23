#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <iostream>
#include <vector>

#include "hist_fast.h"
#include "ConfigParser.h"

using namespace std;
using namespace TMath;

struct Galaxy
{
	double phi;
	double theta;
	double r;       // change to z, allowing new cosmology
	double w;
};

double z2r(const double z)
{
// cosmological parameters
	const double omegaM=0.274;
	const double omegaL=1.-omegaM;
//  conversion to Mpc/h-1
	const double norm=3000.;

	double func=z*(1-omegaM*3*z/4.);
	double r=norm*func;
	return r;
}

// Note: assumes flat space. Would want to replace with
// Dc, Dm to get distance given z and a particular cosmology
inline double dist(const Galaxy& A, const Galaxy& B)
{
	double C = Cos(A.phi)*Sin(A.theta) * Cos(B.phi)*Sin(B.theta) +
             Sin(A.phi)*Sin(A.theta) * Sin(B.phi)*Sin(B.theta) +
             Cos(A.theta) * Cos(B.theta);
	return Sqrt(A.r*A.r + B.r*B.r - 2.*B.r*A.r*C);
}

// Angular map for galaxy positions. Hard-codes rectangular bins
// in theta and phi.
template< typename T > class Map2D
{
	private:
		size_t nxbins_, nybins_;
		double xmin_, xmax_;
		double ymin_, ymax_;
		vector< vector< T > > data_;
	public:
		Map2D(size_t nxbins, double xmin, double xmax, size_t nybins, double ymin, double ymax) : nxbins_(nxbins), nybins_(nybins), xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax), data_(nxbins*nybins)
	{
	}
		int getNumBinsX() const {return nxbins_;}
		int getNumBinsY() const {return nybins_;}
		int getNumBins() const {return nxbins_*nybins_;}
		int getXBin(double x) const {return (x-xmin_)*nxbins_/(xmax_ - xmin_);}
		int getYBin(double y) const {return (y-ymin_)*nybins_/(ymax_ - ymin_);}
		int getBin(double x, double y) const {return getXBin(x) + nxbins_ * getYBin(y);}

		void fill(double x, double y, const T& obj)
		{
			int bin = getBin(x, y);
			if(bin < 0 || bin >= int(data_.size())) return;

			data_[bin].push_back(obj);
		}

		int getBinByBins(int xbin, int ybin) const {return xbin + nxbins_ * ybin;}
		int getXBinByBin(int bin) const {return bin % nxbins_;}
		int getYBinByBin(int bin) const {return bin / nxbins_;}

		size_t getBinEntries(int bin) const
		{
			return data_[bin].size();
		}
		size_t getBinEntries(int xbin, int ybin) const
		{
			return getBinEntries(getBinByBins(xbin, ybin));
		}
		vector<T>& getBinValue(int bin)
		{
			return data_[bin];
		}
		vector<T>& getBinValue(int xbin, int ybin)
		{
			return getBinValue(getBinByBins(xbin, ybin));
		}
		const vector<T>& getBinValue(int bin) const
		{
			return data_[bin];
		}
		const vector<T>& getBinValue(int xbin, int ybin) const
		{
			return getBinValue(getBinByBins(xbin, ybin));
		}
};

// Galaxy 2D angular position and weight, using
// direction cosines of theta and phi
struct Galaxy_ang
{
	double w;
	double cphi;
	double sphi;
	double ctheta;
	double stheta;
};

class Correlations
{
	private:
		double thetamin_;
		double thetamax_;
		double phimin_;
		double phimax_;
		double rmin_;
		double rmax_;
		int thetaregions_;
		int phiregions_;
		int rbins_;
		int abins_;
		int thetabins_;
		int phibins_;
		double smin_;
		double smax_;
		int sbins_;
		string outfile_;

		Map2D<Galaxy>* D;
		Map2D<Galaxy_ang>* R;
		double dw = 0.;
		double dww = 0.;
		double rw = 0.;
		double rww = 0.;

		Hist1D* RR_r = nullptr;	
		Hist1D* RR_alpha = nullptr;
		Hist2D* DR_alpha_r = nullptr;
		Hist1D* DD_cor = nullptr;
		TH1D* htime = nullptr;
		TH1D* hnorm = nullptr;

		void LoadMCTree(string filename)
		{
			Galaxy pos;
			Galaxy_ang pos_ang;

			TFile* tf = TFile::Open(filename.c_str());
			TTree* tr = dynamic_cast<TTree*>(tf->Get("data_pol"));
			tr->SetBranchAddress("position_pol", &pos);
			Hist2D MC_phi_theta(phibins_, phimin_, phimax_, thetabins_, thetamin_, thetamax_);
			rw = 0.;
			rww = 0.;

      // Hack (2016/06/14): convert position theta and r to reflect input:
      // theta is actually declination in input
      // r is actually z in input
      pos.theta = 0.5*Pi() - pos.theta;
      pos.r = z2r(pos.r);

      // Loop over entries in TTree (theta, phi, r, w in each record)
			for(int n = 0 ; n < tr->GetEntries() ; ++n)
			{
				tr->GetEntry(n);
        // Fill 2D distribution.
				MC_phi_theta.fill(pos.phi, pos.theta);

        // Fill radial distribution with weights.
				RR_r->fill(pos.r, pos.w);

        // Track weight sums and sum^2
				rw += pos.w;
				rww += pos.w*pos.w;
			}

      // Create 2D map with the mean positions and weights (random)
			for(int b = 0 ; b < MC_phi_theta.getNumBins() ; ++b)
			{
				if(MC_phi_theta.getBinValue(b) != 0.)
				{
					pos_ang.w = MC_phi_theta.getBinValue(b);
					pos_ang.cphi = Cos(MC_phi_theta.getBinMeanX(b));
					pos_ang.sphi = Sin(MC_phi_theta.getBinMeanX(b));
					pos_ang.ctheta = Cos(MC_phi_theta.getBinMeanY(b));
					pos_ang.stheta = Sin(MC_phi_theta.getBinMeanY(b));
					R->fill(MC_phi_theta.getBinMeanX(b), MC_phi_theta.getBinMeanY(b), pos_ang);
				}
			}
			tf->Close();
		}

		void LoadDATree(string filename)
		{
			Galaxy pos;

			TFile* tf = TFile::Open(filename.c_str());
			TTree* tr = dynamic_cast<TTree*>(tf->Get("data_pol"));
			tr->SetBranchAddress("position_pol", &pos);
			dw = 0.;
			dww = 0;

      // Hack (2016/06/14): convert position theta and r to reflect input:
      // theta is actually declination in input
      // r is actually z in input
      pos.theta = 0.5*Pi() - pos.theta;
      pos.r = z2r(pos.r);

      // Store data galaxy positions, weights, and weight sums
			for(int n = 0 ; n < tr->GetEntries() ; ++n)
			{
				tr->GetEntry(n);
				D->fill(pos.phi, pos.theta, pos);
				dw += pos.w;
				dww += pos.w*pos.w;
			}
			tf->Close();
		}

		void calculateRR(int bina, int binb)
		{
			cout << "Start RR " << bina << " " << binb << endl;

      // Correlations in two different random bins
			if(bina != binb)
			{
				const vector<Galaxy_ang>& va = R->getBinValue(bina);
				const vector<Galaxy_ang>& vb = R->getBinValue(binb);

        // Compute histogram of galaxy-galaxy angular distance,
        // weighted by galaxy weights.
				for(const Galaxy_ang& ga : va)
				{
					for(const Galaxy_ang& gb : vb)
					{
						double dist = ga.cphi*ga.stheta * gb.cphi*gb.stheta;
						dist += ga.sphi*ga.stheta *  gb.sphi*gb.stheta;
						dist += ga.ctheta		  *  gb.ctheta;
						RR_alpha->fill(ACos(dist), ga.w*gb.w);	
					}
				}
			}
      // Correlations within the same bin
			else
			{
				const vector<Galaxy_ang>& va = R->getBinValue(bina);
				for(size_t x = 0 ; x < va.size() ; ++x)
				{   
					const Galaxy_ang& ga = va[x];
					for(size_t y = 0 ; y <= x ; ++y)
					{   
						const Galaxy_ang& gb = va[y];
						double dist = ga.cphi*ga.stheta * gb.cphi*gb.stheta;
						dist += ga.sphi*ga.stheta *  gb.sphi*gb.stheta;
						dist += ga.ctheta		  *  gb.ctheta;
						double f = 1.;
						if(x == y) {f = 0.5;}
						RR_alpha->fill(ACos(dist), f*ga.w*gb.w);	
					}
				}
			}
		}

		void calculateRD(size_t bina, size_t binb)
		{
			cout << "Start RD " << bina << " " << binb << endl;
			const vector<Galaxy>& va = D->getBinValue(bina);
			const vector<Galaxy_ang>& vb = R->getBinValue(binb);
			for(const Galaxy& ga : va)
			{
				double vx = Cos(ga.phi)*Sin(ga.theta);
				double vy = Sin(ga.phi)*Sin(ga.theta);
				double vz = Cos(ga.theta);
				for(const Galaxy_ang& gb : vb)
				{
					double dist = vx *  gb.cphi*gb.stheta;
					dist += vy *  gb.sphi*gb.stheta;
					dist += vz *  gb.ctheta;

          // Note: calculation as a function of r assumes a
          // cosmology. Replace with z?
					DR_alpha_r->fill(ga.r, ACos(dist), ga.w*gb.w);	
				}
			}
		}


		void calculateDD(int bina, int binb)
		{   
			cout << "Start DD " << bina << " " << binb << endl;
			if(bina != binb)
			{
				const vector<Galaxy>& va = D->getBinValue(bina);
				const vector<Galaxy>& vb = D->getBinValue(binb);
				for(const Galaxy& ga : va)
				{
					for(const Galaxy& gb : vb)
					{
            // 1D histogram of DD as a function of 3D separation.
            // Note: distance calculation assumes flat space.
						DD_cor->fill(dist(ga,gb), ga.w*gb.w);
					}
				}
			}
			else
			{
				const vector<Galaxy>& va = D->getBinValue(bina);
				for(size_t x = 0 ; x < va.size() ; ++x)
				{   
					for(size_t y = 0 ; y < x ; ++y)
					{   
            // 1D histogram of DD as a function of 3D separation.
            // Note: distance calculation assumes flat space.
						DD_cor->fill(dist(va[x],va[y]), va[x].w*va[y].w);
					}
				}
			}
		}

    // Parallel job management
		pair<int, int > getjobrange(int job_n, int job_tot, int totev)
		{   
			int nev = 0;
			int startev = 0;
			int perjob = totev/job_tot;
			if(job_n < (totev % job_tot))
			{   
				nev = perjob + 1; 
				startev = (perjob + 1)*job_n;
			}
			else
			{   
				nev = perjob;
				startev = job_n*perjob + (totev % job_tot);
			}
			return pair<int, int>(startev, startev+nev);
		}

	public:
		Correlations(string configfile)
	{
		ConfigParser cfg(configfile);
		outfile_ = cfg.Get<string>("file_out");
		thetamin_ = cfg.Get<double>("theta_min");
		thetamax_ = cfg.Get<double>("theta_max");
		phimin_ = cfg.Get<double>("phi_min");
		phimax_ = cfg.Get<double>("phi_max");
		rmin_ = cfg.Get<double>("r_min");
		rmax_ = cfg.Get<double>("r_max");
		thetaregions_ = cfg.Get<int>("theta_regions");
		phiregions_ = cfg.Get<int>("phi_regions");
		thetabins_ = cfg.Get<int>("theta_bins");
		phibins_ = cfg.Get<int>("phi_bins");
		rbins_ = cfg.Get<int>("r_bins");
		abins_ = cfg.Get<int>("alpha_bins");
		smin_ = cfg.Get<double>("s_min");
		smax_ = cfg.Get<double>("s_max");
		sbins_ = cfg.Get<int>("s_bins");

		R = new Map2D<Galaxy_ang>(phiregions_, phimin_, phimax_, thetaregions_, thetamin_, thetamax_);
		D = new Map2D<Galaxy>(phiregions_, phimin_, phimax_, thetaregions_, thetamin_, thetamax_);
		RR_r = new Hist1D(rbins_, rmin_, rmax_);
		RR_alpha = new Hist1D(abins_, 0, Pi());
		DR_alpha_r = new Hist2D(rbins_, rmin_, rmax_, abins_, 0, Pi());
		DD_cor = new Hist1D(sbins_, smin_, smax_);

		htime = new TH1D("htime", "htime", 10, 0, 10);
		hnorm = new TH1D("hnorm", "hnorm", 3, 0, 3);

		LoadMCTree(cfg.Get<string>("file_random"));
		LoadDATree(cfg.Get<string>("file_data"));
	}

		~Correlations()
		{
			TFile* fout = TFile::Open(outfile_.c_str(), "recreate");
			RR_r->writeTH1D("RR_r");
			RR_alpha->writeTH1D("RR_alpha");
			DR_alpha_r->writeTH2D("DR_alpha_r");
			DD_cor->writeTH1D("DD_cor");
			htime->Write("htime");
			hnorm->Write("hnorm");
			fout->Write();
			fout->Close();
		}

		void Calculate(int job_n, int job_tot)
		{
			stringstream ss;
			ss << outfile_ << "_" << job_n << "_" << job_tot << ".root";
			outfile_ = ss.str();
			if(job_n == 0)
			{
				hnorm->SetBinContent(1, (rw*rw-rww)*0.5);
				hnorm->SetBinContent(2, rw*dw);
				hnorm->SetBinContent(3, (dw*dw-dww)*0.5);
				cout << "Normalization RR = " <<  hnorm->GetBinContent(1) << ", RD = " << hnorm->GetBinContent(2) << ", DD = " << hnorm->GetBinContent(3) << endl; 
			}
			int av_jobs = R->getNumBins()*5;
			cout << "Job number: " << job_n << ", total number of jobs: " << job_tot << ", available jobs: " << av_jobs << endl;
			pair<int, int > range(getjobrange(job_n, job_tot, av_jobs));
			cout << range.first << " " << range.second << endl;
			for(int i = range.first ; i < range.second ; ++i)
			{
				int bina = i / 5;
				int nb = i % 5;
				int xb = R->getXBinByBin(bina);
				int yb = R->getYBinByBin(bina);
				if(nb == 1){ xb++;}	
				else if(nb == 2){ xb--; yb++;}	
				else if(nb == 3){ yb++;}	
				else if(nb == 4){ xb++; yb++;}
				int binb = R->getBinByBins(xb, yb);
				int start;
				if(bina > 0 && binb > 0 && binb < R->getNumBins() && bina < R->getNumBins())
				{
					start = time(nullptr);	
					calculateDD(bina, binb);
					htime->Fill(0.5, time(nullptr) - start);
					start = time(nullptr);	
					calculateRD(bina, binb);
					if(bina != binb) {calculateRD(binb, bina);};
					htime->Fill(1.5, time(nullptr) - start);
					start = time(nullptr);	
					calculateRR(bina, binb);
					htime->Fill(2.5, time(nullptr) - start);
				}
			}
		}

};


int main(int argc, char** argv)
{
    int job_n = atoi(argv[1]);
    int job_tot = atoi(argv[2]);
    string configfile(argv[3]);

    Correlations cor(configfile);
    cor.Calculate(job_n, job_tot);
}
