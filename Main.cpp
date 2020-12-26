#include "Particle.h"

#include "TH1F.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"

int main()
{
    auto start = std::chrono::system_clock::now();



    int nBins = 1e3;
    int nGen = 1e5;
    int N = 120;
    int NParticle = 100;
    double phi;
    double theta;
    double impulse;

    char p_plus[] = "Pioni+";
    char p_minus[] = "Pioni-";
    char k_plus[] = "Kaoni+";
    char k_minus[] = "Kaoni-";
    char pr_plus[] = "Protoni+";
    char pr_minus[] = "Protoni-";
    char Kstar[] = "K*";

    double constexpr pMass = 0.13957;
    double constexpr kMass = 0.49367;
    double constexpr prMass = 0.93827;
    double constexpr KstarMass = 0.89166;

    Particle::AddParticle(p_plus, pMass, 1);  //ID = 1
    Particle::AddParticle(p_minus, pMass, -1);  // ID = 2
    Particle::AddParticle(k_plus, kMass, 1);  // ID = 3
    Particle::AddParticle(k_minus, kMass, -1);  // ID = 4
    Particle::AddParticle(pr_plus, prMass, 1);  // ID = 5
    Particle::AddParticle(pr_minus, prMass, -1);  // ID = 6
    Particle::AddParticle(Kstar, KstarMass, 0, 0.050);  // ID = 7
    std::cout << "Particle added successfully" << '\n';


    TFile *file = new TFile("data.root", "recreate");
    std::cout << "File created successfully " << '\n';


    TH1F *h1 = new TH1F("h1", "Particle's type", 7, 0, 7);
    TH1F *h2 = new TH1F("h2", "Azimuthal angle distribution", 500, 0, 2 * M_PI);
    TH1F *h3 = new TH1F("h3", "Polar angle distribution", 250, 0, M_PI);
    TH1F *h4 = new TH1F("h4", "Impulse", nBins, 0, 5.5);
    TH1F *h5 = new TH1F("h5", "Transverse impulse", nBins, 0, 4); //senza componente z
    TH1F *h6 = new TH1F("h6", "Energy", nBins, 0, 5);
    TH1F *h7 = new TH1F("h7", "Invariant masses", nBins, 0, 5);
    TH1F *h8 = new TH1F("h8", "Invariant m. different charge", nBins, 0, 5);
    TH1F *h9 = new TH1F("h9", "Invariant m. same charge", nBins, 0, 5);
    TH1F *h10 = new TH1F("h10", "P+/K- , P-/K+", nBins, 0, 5);
    TH1F *h11 = new TH1F("h11", "P+/K+, P-/K-", nBins, 0, 5);
    TH1F *h12 = new TH1F("h12", "K* decays", nBins, 6, 1.3);
    std::cout << "Histograms created successfully" << '\n';

    h7->Sumw2();
    h8->Sumw2();
    h9->Sumw2();
    h10->Sumw2();
    h11->Sumw2();
    h12->Sumw2();


    gRandom -> SetSeed();


    std::cout << "Beginning " << nGen << " generations, please wait" << '\n'; 
    for(int i = 0; i!= nGen; ++i)
    {
        int particleTail = 0;
        Particle particle[N];
        for(int j = 0; j != NParticle; ++j)
        {
            phi = gRandom -> Uniform(2 * M_PI);
            theta = gRandom -> Uniform(M_PI);
            h2 -> Fill(phi);
            h3 -> Fill(theta);

            impulse = gRandom -> Exp(1);
            double px = impulse * sin(theta) * cos(phi);
            double py = impulse * sin(theta) * sin(phi);
            double pz = impulse * cos(theta);
            particle[j].SetImpulse(px, py, pz);
            h4 -> Fill(impulse);
            h5 -> Fill(sqrt(px * px + py * py));

            double index = gRandom -> Uniform(1);
            if(index < 0.4)     //80%   pioni
            {
                particle[j].SetParticleName(p_plus);
                h1 -> Fill(0);
            }
            else if(index < 0.8)     //80%
            {
                particle[j].SetParticleName(p_minus);
                h1 -> Fill(1);
            }
            else if(index < 0.85)    //10%   kaoni
            {
                particle[j].SetParticleName(k_plus);
                h1 -> Fill(2);
            }
            else if(index < 0.9)    //10%
            {
                particle[j].SetParticleName(k_minus);
                h1 -> Fill(3);
            }
            else if(index < 0.945)    //9%   protoni
            {
                particle[j].SetParticleName(pr_plus);
                h1 -> Fill(4);
            }
            else if(index < 0.99)    //9%
            {
                particle[j].SetParticleName(pr_minus);
                h1 -> Fill(5);
            }
            else if(index < 0.995)    //1%   K risonanze
            {
                particle[j].SetParticleName(Kstar);
                h1 -> Fill(6);
                particle[NParticle + particleTail].SetParticleName(p_plus);
                particle[NParticle + particleTail + 1].SetParticleName(k_minus);
                particle[j].Decay2body(particle[NParticle + particleTail + 1], particle[NParticle + particleTail]);
                h12 -> Fill(particle[NParticle + particleTail + 1].GetInvMass(particle[NParticle + particleTail]));
                ++particleTail;
                ++particleTail;
            }
            else    //1%
            {
                particle[j].SetParticleName(Kstar);
                h1 -> Fill(6);
                particle[NParticle + particleTail].SetParticleName(p_minus);
                particle[NParticle + particleTail + 1].SetParticleName(k_plus);
                particle[j].Decay2body(particle[NParticle + particleTail + 1], particle[NParticle + particleTail]);
                h12 -> Fill(particle[NParticle + particleTail + 1].GetInvMass(particle[NParticle + particleTail]));
                ++particleTail;
                ++particleTail;
            }
            h6 -> Fill(particle[j].GetParticleEnergy());
        }
        
        for(int g = 0; g != (NParticle + particleTail - 1); ++g)
        {
            for(int h = (g + 1); h != (NParticle + particleTail); ++h)
            {
                double gID = particle[g].GetIndex();
                double hID = particle[h].GetIndex();
                double MassInvariant = particle[g].GetInvMass(particle[h]);
                bool Couple = ((gID == 0 || gID == 1) && (hID == 2 || hID == 3)); // caso P/K 
                h7 -> Fill(MassInvariant);
                
                if((gID == 0) || (gID == 2) || (gID == 4))
                {
                    if((hID == 0) || (hID == 2) || (hID == 4))
                    {
                        h9 -> Fill(MassInvariant);
                        if(Couple)
                        {
                            h11 -> Fill(MassInvariant);
                        }
                    }
                    else if((hID == 1) || (hID == 3) || (hID == 5))
                    {
                        h8 -> Fill(MassInvariant);
                        if(Couple)
                        {
                            h10 -> Fill(MassInvariant);
                        }
                    }
                }
                else if((gID == 1) || (gID == 3) || (gID == 5))
                {
                    if((hID == 0) || (hID == 2) || (hID == 4))
                    {
                        h8 -> Fill(MassInvariant);
                        if(Couple)
                        {
                            h10 -> Fill(MassInvariant);
                        }
                    }
                    else if((hID == 1) || (hID == 3) || (hID == 5))
                    {
                        h9 -> Fill(MassInvariant);
                        if(Couple)
                        {
                            h11 -> Fill(MassInvariant);
                        }
                    }
                }
            }
        }
    }
    std::cout << "Generation of " << nGen << " events ended successfully" << '\n';
    
    file -> Write();
    std::cout << "File transcripted successfully " << '\n';
    file -> Close();

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> diff = end-start;
    std::cout << "Time required: " << diff.count() << " s\n";


    return 0;
}
