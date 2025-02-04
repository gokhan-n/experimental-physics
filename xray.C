{
    double d = 5.64*1e-10;
    double rad = 180/(TMath::Pi());
    double slope,intercept;
    double slopeW,interceptW,slopeW_error,interceptW_error; 
    float eV = 1.60217663e-19;
    double c = 299792458.;
    
    std::vector<std::string> filenames = {
        "calib35kv.txt",
        "15kv.txt",
        "18kv.txt",
        "21kv.txt",
        "24kv.txt",
        "27kv.txt",
        "30kv.txt",
        

     };

        std::vector<std::string> titles = {
        "calib35v",
        "15kv",
        "18kv",
        "21kv",
        "24kv",
        "27kv",
        "30kv",

        
    };

       std::vector<std::vector<double>> peakRanges = {
        {7.38,7.69},{13.2,13.51},{14.85,15.29},{19.85,20.25},{22.3,22.86}  


    };

    std::vector<std::vector<double>> linefitRanges = {
        {0.172*1e-9,0.187*1e-9},
        {0.142*1e-9,0.156*1e-9},
        {0.125*1e-9,0.14*1e-9},
        {0.118*1e-9,0.138*1e-9},
        {0.10*1e-9,0.115*1e-9} ,
        {0.095*1e-9,0.113*1e-9},
    };


    float peakValues[5],sigmaValues[5]; 
    float theoreticalPeaks[5]= {7.76,13.93,15.69,19.73,23.9};

    float mumins[6],musigmas[6],fmaxs[6],fsigmas[6]; 

    float energies[6] = {eV*15000,eV*18000,eV*21000,eV*24000,eV*27000,eV*30000};



    for(int j =0;j<filenames.size();j++){
 
      

    if (j==0){
    TGraphErrors *gr = new TGraphErrors(filenames[j].c_str());
    gr->Draw("A*");
    gr->GetXaxis()->SetTitle("Angle (Degree)");
    gr->GetYaxis()->SetTitle("Accumulation Rate(Impulse/Second)");
    gr->SetTitle("Accumulation Rate vs Angle");

    for (int i = 0; i < 5; ++i) {
            TF1 *gauss = new TF1(Form("gauss_%d_%d", j, i), "gaus", peakRanges[i][0], peakRanges[i][1]);
            gr->Fit(gauss, "R");


            gauss->Draw("same");
            double peak = gauss->GetParameter(1);
            double sigma = gauss->GetParameter(2);
            cout << "the "<< i <<  "'th  peak is" << peak <<"and its sigma is="<< sigma << "\n";
            peakValues[i]=peak;            
            sigmaValues[i]=sigma;
            

    

        }
        TCanvas *c = new TCanvas();


        TGraphErrors *calib = new TGraphErrors(5,peakValues,theoreticalPeaks,sigmaValues,0);
        calib->Draw("A*");
        calib->GetXaxis()->SetTitle("Experimental Theta Values(Degrees)");
        calib->GetYaxis()->SetTitle("Nominal Thteta Values(Degrees)");
        calib->SetTitle("Experimental vs Nominal Theta Values");

        TF1 *linefit = new TF1("linefit","[0]*x+[1]");
        linefit->SetParNames("Slope", "Intercept");
        linefit->SetParameters(1,0); 
        calib->Fit(linefit);
        gStyle->SetOptFit(1);

        slope = linefit->GetParameter(0);
        intercept = linefit->GetParameter(1);

        TCanvas *c2 = new TCanvas();


    }
    
    else if (j != 0) {
    int ndata;

    std::ifstream inputFile(filenames[j].c_str());
    std::string line;
    int linecount = 0;

    // First, count the number of lines in the file
    while (std::getline(inputFile, line)) {
        linecount++;
    }

    ndata = linecount;

    // Now, rewind the file stream to read data
    inputFile.clear();
    inputFile.seekg(0, std::ios::beg);

    // Create vectors to store the data
    std::vector<float> angles, impulses, angleErrors, impulseErrors;

    float angle, impulse, angleError, impulseError;
    while (std::getline(inputFile, line)) {
        std::istringstream iss(line);
        iss >> angle >> impulse >> angleError >> impulseError;
        
        // Push the values into the vectors
        angles.push_back(angle);
        impulses.push_back(impulse);
        angleErrors.push_back(angleError);
        impulseErrors.push_back(impulseError);
    }

    ndata = angles.size();
    std::vector<float> wavelengths(ndata), wavelengthErrors(ndata);

    for (int i = 0; i < ndata; ++i) {
        wavelengths[i] = 2 * d * sin((angles[i] * slope + intercept) / rad);
        wavelengthErrors[i] = 2 * d * cos((angles[i] * slope + intercept) / rad) * (0.1 / rad);
    }

    TGraphErrors *gr = new TGraphErrors(ndata, wavelengths.data(), impulses.data(), wavelengthErrors.data(), impulseErrors.data());
    gr->Draw("A*");
    gr->GetXaxis()->SetTitle("Wavelength(m)");
    gr->GetYaxis()->SetTitle("Accumulation Rate(Impulse/Second)");
    gr->SetTitle(Form("Accumulation Rate vs Wavelength - %s", titles[j].c_str()));

    TF1 *linefitW = new TF1("linefit","[0]*x+[1]");
    linefitW->SetParNames("Slope", "Intercept");
    linefitW->SetParameters(1.15*1e12,-200); 
    linefitW->SetRange(linefitRanges[j - 1][0], linefitRanges[j - 1][1]); // Set the fitting range


    gr->Fit(linefitW,"R");
    slopeW = linefitW->GetParameter(0);
    interceptW = linefitW->GetParameter(1);

    slopeW_error = linefitW->GetParError(0);
    interceptW_error = linefitW->GetParError(1);

    cout << "the slope is  "<< slopeW<<  "'+- "<< slopeW_error <<  "'the intercept is "
     <<interceptW<<  "'+- "<< interceptW_error<<  "\n";

     mumins[j-1] = -interceptW / slopeW;
    musigmas[j-1] = sqrt(
        pow(interceptW_error / slopeW, 2) +
        pow(interceptW * slopeW_error / (slopeW * slopeW), 2)
    );

    fmaxs[j-1] = c/mumins[j-1];
    fsigmas[j-1] = (c * musigmas[j-1]) / (mumins[j-1] * mumins[j-1]);

    cout << "mu mins are  "<< mumins[j-1] <<  " +- " <<musigmas[j-1]<< "\n";
    cout << "f max are  "<< fmaxs[j-1] <<  " +- " <<fsigmas[j-1]<< "\n";




    TCanvas *c3 = new TCanvas();
}
    }

TGraphErrors *planck = new TGraphErrors(6,fmaxs,energies,0,0);
    planck->Draw("A*");
    planck->GetXaxis()->SetTitle("Frequency(Hz)");
    planck->GetYaxis()->SetTitle("Energy(J)");
    planck->SetTitle("Energy vs Frequency");

     TF1 *linePlanck = new TF1("linefit","[0]*x+[1]");
        linePlanck->SetParNames("Slope", "Intercept");
        linePlanck->SetParameters(6*1e-34,1*1e-15); 
        planck->Fit(linePlanck);

        double Planck = linePlanck->GetParameter(0);
        double PlanckError = linePlanck->GetParError(0);

        cout << "planck constant is found out to be  "<< Planck<< " +- " <<PlanckError<<"\n";


}