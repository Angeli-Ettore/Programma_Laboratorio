void Analisis()
{
    TFile *file = new TFile("data.root");
    
    gStyle -> SetOptFit(110);
    gStyle -> SetOptStat(210);

    TCanvas *c1 = new TCanvas();
    TCanvas *c2 = new TCanvas();
    TCanvas *c3 = new TCanvas();
    TCanvas *c4 = new TCanvas();
    c1 -> Divide(2, 2);
    c2 -> Divide(3, 2);
    c3 -> Divide(2);
    c4 -> Divide(3);

    TH1F *k[12];
    k[0] = (TH1F*)file -> Get("h1");
    k[1] = (TH1F*)file -> Get("h2");
    k[2] = (TH1F*)file -> Get("h3");
    k[3] = (TH1F*)file -> Get("h4");
    k[4] = (TH1F*)file -> Get("h5");
    k[5] = (TH1F*)file -> Get("h6");
    k[6] = (TH1F*)file -> Get("h7");
    k[7] = (TH1F*)file -> Get("h8");
    k[8] = (TH1F*)file -> Get("h9");
    k[9] = (TH1F*)file -> Get("h10");
    k[10] = (TH1F*)file -> Get("h11");
    k[11] = (TH1F*)file -> Get("h12");

    
    //K COSMETIC
    k[0] -> GetXaxis() -> SetTitle("p+/p-           k+/k-           pr+/pr-           k*");
    k[1] -> GetXaxis() -> SetTitle("Azimuthal Angle (rad)");
    k[2] -> GetXaxis() -> SetTitle("Polar Angle (rad)");
    k[3] -> GetXaxis() -> SetTitle("Momentum (Gev)");
    k[4] -> GetXaxis() -> SetTitle("Transverse momentum (Gev)");
    k[5] -> GetXaxis() -> SetTitle("Energy");


    for(int i = 0; i != 12; ++i)
    {
        k[i] -> SetFillColor(35);
        k[i] -> SetMarkerStyle(kFullDotLarge);
        k[i] -> SetMarkerColor(40);
        k[i] -> SetMarkerSize(0.20);
    }

    for(int j = 0; j != 4; ++j)
    {
        k[j] -> GetXaxis() -> SetLabelFont(82);
        k[j] -> GetXaxis() -> SetTitleFont(82);
        k[j] -> GetXaxis() -> SetTitleOffset(1);
        k[j] -> GetYaxis() -> SetTitle("Occurrences");
        k[j] -> GetYaxis() -> SetLabelFont(82);
        k[j] -> GetYaxis() -> SetTitleFont(82);
        k[j] -> GetYaxis() -> SetTitleOffset(1.5);
    }

    for(int g = 1; g != 3; ++g)
    {
        c3 -> cd(g);
        k[g + 3] -> GetXaxis() -> SetLabelFont(82);
        k[g + 3] -> GetXaxis() -> SetTitleFont(82);
        k[g + 3] -> GetXaxis() -> SetTitleOffset(1);
        k[g + 3] -> GetYaxis() -> SetTitle("Occurrences");
        k[g + 3] -> GetYaxis() -> SetLabelFont(82);
        k[g + 3] -> GetYaxis() -> SetTitleFont(82);
        k[g + 3] -> GetYaxis() -> SetTitleOffset(1.5);
        k[g + 3] -> DrawCopy();
    }
    
    for(int f = 1; f != 7; ++f)
    {
        c2 -> cd(f);
        k[f + 5] -> GetXaxis() -> SetTitle("Mass Invariant (Gev/c^2)");
        k[f + 5] -> GetXaxis() -> SetLabelFont(82);
        k[f + 5] -> GetXaxis() -> SetTitleFont(82);
        k[f + 5] -> GetXaxis() -> SetTitleOffset(1);
        k[f + 5] -> GetYaxis() -> SetTitle("Occurrences");
        k[f + 5] -> GetYaxis() -> SetLabelFont(82);
        k[f + 5] -> GetYaxis() -> SetTitleFont(82);
        k[f + 5] -> GetYaxis() -> SetTitleOffset(1);
        k[f + 5] -> DrawCopy();
    }



//Checking proportions

    c1 -> cd(1);
    k[0] -> DrawCopy();

    double total = k[0] -> GetEntries();
    double p_plus = k[0] -> GetBinContent(1);
    double p_minus = k[0] -> GetBinContent(2);
    double k_plus = k[0] -> GetBinContent(3);
    double k_minus = k[0] -> GetBinContent(4);
    double pr_plus = k[0] -> GetBinContent(5);
    double pr_minus = k[0] -> GetBinContent(6);
    double Kstar = k[0] -> GetBinContent(7);
    std::cout << "Number of Pions +: " << (p_plus / total) << " ± " << (k[0] -> GetBinError(1) / total) << " (expectation: 40%)"<< '\n';
    std::cout << "Number of Pions -: " << (p_minus / total) << " ± " << (k[0] -> GetBinError(2) / total) << " (expectation: 40%)"<< '\n';
    std::cout << "Number of Kaons +: " << (k_plus / total) << " ± " << (k[0] -> GetBinError(3) / total) << " (expectation: 5%)"<< '\n';
    std::cout << "Number of Kaons -: " << (k_minus / total) << " ± " << (k[0] -> GetBinError(4) / total) << " (expectation: 5%)"<< '\n';
    std::cout << "Number of Protons +: " << (pr_plus / total) << " ± " << (k[0] -> GetBinError(5) / total) << " (expectation: 4.5%)"<< '\n';
    std::cout << "Number of Protons -: " << (pr_minus / total) << " ± " << (k[0] -> GetBinError(6) / total) << " (expectation: 4.5%)"<< '\n';
    std::cout << "Number of K* : " << (Kstar / total) << " ± " << (k[0] -> GetBinError(7) / total) << " (expectation: 1%)"<< '\n' << '\n';



//Angle fit
    c1 -> cd(2);
    k[1] -> Fit("pol0", "Q");
    k[1] -> DrawCopy();

    c1 -> cd(3);
    k[2] -> Fit("pol0", "Q");
    k[2] -> DrawCopy();

    TF1 *fit1 = k[1] -> GetFunction("pol0");
    TF1 *fit2 = k[2] -> GetFunction("pol0");

    std::cout << "Fitting Angle's distribution: " << '\n';
    double Chi1 = fit1 -> GetChisquare();
    double Dof1 = fit1 -> GetNDF();
    std::cout << "Azimutal Angle: " << '\n';
    std::cout << "Chi Square: " << Chi1 << '\n';
    std::cout << "Degrees of freedom: " << Dof1 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi1 / Dof1) << '\n'<< '\n';


    double Chi2 = fit2 -> GetChisquare();
    double Dof2 = fit2 -> GetNDF();
    std::cout << "Polar Angle: " << '\n';
    std::cout << "Chi Square: " << Chi2 << '\n';
    std::cout << "Degrees of freedom: " << Dof2 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi2 / Dof2) << '\n'<< '\n';



//Impulse fit
    c1 -> cd(4);
    k[3] -> Fit("expo", "Q");
    k[3] -> DrawCopy();
    

    TF1 *fit3 = k[3] -> GetFunction("expo");

    std::cout << "Fitting Impulse's distribution: " << '\n';
    double Chi3 = fit3 -> GetChisquare();
    double Dof3 = fit3 -> GetNDF();
    double Mean3 = k[3] -> GetMean();
    double MeanError3 = k[3] -> GetMeanError();
    std::cout << "Distribution's mean: " << Mean3 << " ± " << MeanError3 << " (expectation: 1)" << '\n';
    std::cout << "Chi Square: " << Chi3 << '\n';
    std::cout << "Degrees of freedom: " << Dof3 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi3 / Dof3) << '\n' << '\n';



//K* fit
    gStyle -> SetOptStat(210);
    c4 -> cd(1);
    k[11] -> Fit("gaus", "Q");
    k[11] -> DrawCopy();

    TF1 *fit4 = k[11] -> GetFunction("gaus");

    std::cout << "Fitting K* Decay distribution: " << '\n';
    double Chi4 = fit4 -> GetChisquare();
    double Dof4 = fit4 -> GetNDF();
    double Mean4 = fit4 -> GetParameter(1);
    double MeanError4 = fit4 -> GetParError(1);
    double Sigma4 = fit4 -> GetParameter(2);
    double SigmaError4 = fit4 -> GetParError(2);
    std::cout << "Distribution's mean: " << Mean4 << " ± " << MeanError4 << '\n';
    std::cout << "Distribution's RMS: " << Sigma4 << " ± " << SigmaError4 << '\n';
    std::cout << "Chi Square: " << Chi4 << '\n';
    std::cout << "Degrees of freedom: " << Dof4 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi4 / Dof4) << '\n' << '\n';

//Subtracting PK concordant to PK discordant
    TH1F kSubtractionPK1 = (*k[9] - *k[10]);
    kSubtractionPK1.SetTitle("Subtraction between discordant PK and concordant PK");
    kSubtractionPK1.SetAxisRange(0.0, 2.0);
    kSubtractionPK1.GetYaxis() -> SetTitleFont(82);
    kSubtractionPK1.GetYaxis() -> SetTitleOffset(1);
    
    TF1 *fit5 = new TF1("Fit5", "gaus", 0.6, 1.3);
    TCanvas *c5 = new TCanvas();

    c5 -> Divide(2, 2);
    c5 -> cd(1);
    k[9] -> DrawCopy();
    c5 -> cd(2);
    k[10] -> DrawCopy();
    c5 -> cd(3);
    kSubtractionPK1.Fit("Fit5", "Q");
    kSubtractionPK1.DrawCopy();
    c5 -> cd(4);
    k[11] -> Fit("gaus", "Q");
    k[11] -> DrawCopy();
    c4 -> cd(2);
    kSubtractionPK1.DrawCopy();


    std::cout << "Fitting the subtraction between P-K with discordant and concordant charge" << '\n';
    double Chi5 = fit5 -> GetChisquare();
    double Dof5 = fit5 -> GetNDF();
    double Mean5 = fit5 -> GetParameter(1);
    double MeanError5 = fit5 -> GetParError(1);
    double Sigma5 = fit5 -> GetParameter(2);
    double SigmaError5 = fit5 -> GetParError(2);
    std::cout << "Distribution's mean: " << Mean5 << " ± " << MeanError5 << '\n';
    std::cout << "Distribution's mean from simulation: " << Mean4 << " ± " << MeanError4 << '\n';
    std::cout << "Distribution's RMS: " << Sigma5 << " ± " << SigmaError5 << '\n';
    std::cout << "Distribution's RMS from simulation: " << Sigma4 << " ± " << SigmaError4 << '\n';
    std::cout << "Chi Square: " << Chi5 << '\n';
    std::cout << "Degrees of freedom: " << Dof5 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi5 / Dof5) << '\n' << '\n';


//Subtracting concordant and discordant particle
    TH1F kSubtractionPK2 = (*k[7] - *k[8]);
    kSubtractionPK2.SetTitle("Subtraction between discordant and concordant particles");
    kSubtractionPK2.SetAxisRange(0.0, 2.0);
    kSubtractionPK2.GetYaxis() -> SetTitleFont(82);
    kSubtractionPK2.GetYaxis() -> SetTitleOffset(1);

    TF1 *fit6 = new TF1("Fit6", "gaus", 0.6, 1.3);

    TCanvas *c6 = new TCanvas();
    c6 -> Divide(2, 2);
    c6 -> cd(1);
    k[7] -> DrawCopy();
    c6 -> cd(2);
    k[8] -> DrawCopy();
    c6 -> cd(3);
    kSubtractionPK2.Fit("Fit6", "Q");
    kSubtractionPK2.DrawCopy();
    c6 -> cd(4);
    k[11] -> Fit("gaus", "Q");
    k[11] -> DrawCopy();
    c4 -> cd(3);
    kSubtractionPK2.DrawCopy();


    std::cout << "Fitting the subtraction between all particles with discordant and concordant charge" << '\n';
    double Chi6 = fit6 -> GetChisquare();
    double Dof6 = fit6 -> GetNDF();
    double Mean6 = fit6 -> GetParameter(1);
    double MeanError6 = fit6 -> GetParError(1);
    double Sigma6 = fit6 -> GetParError(2);
    double SigmaError6 = fit6 -> GetParError(2);
    std::cout << "Distribution's mean: " << Mean6 << " ± " << MeanError6 << '\n';
    std::cout << "Distribution's mean from simulation: " << Mean4 << " ± " << MeanError4 << '\n';
    std::cout << "Distribution's RMS: " << Sigma6 << " ± " << SigmaError6 << '\n';
    std::cout << "Distribution's RMS from simulation: " << Sigma4 << " ± " << SigmaError4 << '\n';
    std::cout << "Chi Square: " << Chi6 << '\n';
    std::cout << "Degrees of freedom: " << Dof6 << '\n';
    std::cout << "Reduced Chi Square: " << (Chi6 / Dof6) << '\n' << '\n';



//Printig all the canvas
    c1 -> Print("c1.pdf", "pdf");
    c2 -> Print("c2.pdf", "pdf");
    c3 -> Print("c3.pdf", "pdf");
    c4 -> Print("c4.pdf", "pdf");
}
