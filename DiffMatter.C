 {

 const int ndata = 10;
 const float L=14.0 * 1e-2;
 const float R= 4.3 * 1e-2;
 const float sigmaR = 0.1*1e-2;
 const float sigmaL = 0.3*1e-2;
 const float sigmaArc = 0.01*1e-2;
 const float sigmaV = 100;
double h = 6.62607015 * 1e-34; // exact
double q = 1.602176634 * 1e-19; // exact
double m = 9.1093837015 * 1e-31;


 
 vector<double> voltages = {     
   2800,3000,3200,3400,3600,3800,4100,4300,4500,4700
    };
 
 std::vector<std::vector<double>> mFirst = {
    {0.015, 0.017, 0.014, 0.016},
    {0.014, 0.016, 0.013, 0.015},
    {0.014, 0.016, 0.014, 0.016},
    {0.013, 0.016, 0.013, 0.015},
    {0.013, 0.015, 0.012, 0.014},
    {0.013, 0.014, 0.012, 0.014},
    {0.012, 0.014, 0.013, 0.014},
    {0.012, 0.014, 0.011, 0.013},
    {0.011, 0.013, 0.011, 0.013},
    {0.011, 0.014, 0.011, 0.013}
};

std::vector<std::vector<double>> mSecond = {
    {0.028, 0.030, 0.026, 0.029},
    {0.026, 0.028, 0.025, 0.027},
    {0.025, 0.027, 0.025, 0.027},
    {0.024, 0.027, 0.024, 0.026},
    {0.023, 0.026, 0.023, 0.026},
    {0.023, 0.025, 0.023, 0.026},
    {0.022, 0.024, 0.022, 0.024},
    {0.021, 0.023, 0.021, 0.023},
    {0.020, 0.022, 0.020, 0.022},
    {0.020, 0.023, 0.019, 0.022}
};



float arcsFirst[10],arcsSecond[10], sigmaThetaFirst[10],sigmaThetaSecond[10],
thetasFirst[10],thetasSecond[10],xAxisFirst[10],yAxis[10],xAxisSecond[10],
sigmaxFirst[10],sigmaxSecond[10],sigmaY[10];

for (int i = 0; i < 10; ++i) {
    arcsFirst[i]=((mFirst[i][0]+mFirst[i][1]+mFirst[i][2]+mFirst[i][3])/8);
    arcsSecond[i]=((mSecond[i][0]+mSecond[i][1]+mSecond[i][2]+mSecond[i][3])/8);

    thetasFirst[i]=(atan((sin(arcsFirst[i]/R)*R)/(L-R+R*cos(arcsFirst[i]/R)))/2);
    thetasSecond[i]=(atan((sin(arcsSecond[i]/R)*R)/(L-R+R*cos(arcsSecond[i]/R)))/2);

}

for (int k = 0; k < 10; ++k) {


    double partialR1 = (-(arcsFirst[k] * cos(arcsFirst[k] / R)) / (R * (L - R + R * cos(arcsFirst[k] / R))) 
            + sin(arcsFirst[k] / R) / (L - R + R * cos(arcsFirst[k] / R)) 
            - (R * sin(arcsFirst[k] / R) * (-1 + cos(arcsFirst[k] / R) + (arcsFirst[k] * sin(arcsFirst[k] / R)) / R)) 
            / pow((L - R + R * cos(arcsFirst[k] / R)), 2)) 
           / (2 * (1 + (pow(R, 2) * pow(sin(arcsFirst[k] / R), 2)) / pow((L - R + R * cos(arcsFirst[k] / R)), 2)));



    double partialL1 = -(R * sin(arcsFirst[k] / R)) / 
           (2 * pow(L - R + R * cos(arcsFirst[k] / R), 2) * 
           (1 + (pow(R, 2) * pow(sin(arcsFirst[k] / R), 2)) / pow(L - R + R * cos(arcsFirst[k] / R), 2)));



    double partialArc1 = ((L - R) * cos(arcsFirst[k] / R) + R) / 
           (2 * ((L - R) * (L - R) + 
           2 * (L - R) * R * cos(arcsFirst[k] / R) + 
           R * R * cos(arcsFirst[k] / R) * cos(arcsFirst[k] / R) + 
           R * R * sin(arcsFirst[k] / R) * sin(arcsFirst[k] / R)));




    double partialR2 = (-(arcsSecond[k] * cos(arcsSecond[k] / R)) / (R * (L - R + R * cos(arcsSecond[k] / R))) 
            + sin(arcsFirst[k] / R) / (L - R + R * cos(arcsSecond[k] / R)) 
            - (R * sin(arcsSecond[k] / R) * (-1 + cos(arcsSecond[k] / R) + (arcsSecond[k] * sin(arcsSecond[k] / R)) / R)) 
            / pow((L - R + R * cos(arcsSecond[k] / R)), 2)) 
           / (2 * (1 + (pow(R, 2) * pow(sin(arcsSecond[k] / R), 2)) / pow((L - R + R * cos(arcsSecond[k] / R)), 2)));

    double partialL2 = -(R * sin(arcsSecond[k] / R)) / 
           (2 * pow(L - R + R * cos(arcsSecond[k] / R), 2) * 
           (1 + (pow(R, 2) * pow(sin(arcsSecond[k] / R), 2)) / pow(L - R + R * cos(arcsSecond[k] / R), 2)));

    double partialArc2 = ((L - R) * cos(arcsSecond[k] / R) + R) / 
           (2 * ((L - R) * (L - R) + 
           2 * (L - R) * R * cos(arcsSecond[k] / R) + 
           R * R * cos(arcsSecond[k] / R) * cos(arcsSecond[k] / R) + 
           R * R * sin(arcsSecond[k] / R) * sin(arcsSecond[k] / R)));



    double sigmaTheta1 = sqrt((partialR1*sigmaR*partialR1*sigmaR)+(partialL1*sigmaL*partialL1*sigmaL)+
    (partialArc1*sigmaArc*partialArc1*sigmaArc));
    sigmaThetaFirst[k]=(sigmaTheta1);

    double sigmaTheta2 = sqrt((partialR2*sigmaR*partialR2*sigmaR)+(partialL2*sigmaL*partialL2*sigmaL)+
    (partialArc2*sigmaArc*partialArc2*sigmaArc));
    sigmaThetaSecond[k]=(sigmaTheta2);

}



for (int n = 0; n < 10; ++n) {
  sigmaxFirst[n]=cos(thetasFirst[n])*sigmaThetaFirst[n];
  sigmaxSecond[n] = (cos(thetasSecond[n])*sigmaThetaSecond[n]);
  sigmaY[n] = ((1/pow(voltages[n],1.5))*sigmaV);

  xAxisFirst[n] = sin(thetasFirst[n]);
  xAxisSecond[n]= sin(thetasSecond[n]);

  yAxis[n] = 1/(pow(voltages[n],0.5));

}

TGraphErrors *mygraphFirst = new TGraphErrors(ndata,xAxisFirst,yAxis,sigmaxFirst,sigmaY);
TGraphErrors *mygraphSecond = new TGraphErrors(ndata,xAxisSecond,yAxis,sigmaxSecond,sigmaY);

mygraphFirst->Draw("A*");
TCanvas *c1 = new TCanvas();
TF1 *ffirst = new TF1("ffirst","[0]*x+[1]",0.021,0.029);
mygraphFirst->Fit(ffirst);
double slopeFirst = ffirst->GetParameter(0);
double slope_error_first = ffirst->GetParError(0);
double d_first = slopeFirst * h / sqrt(8 * q * m);
double d_error_first = slope_error_first * h / sqrt(8 *
q * m);
 mygraphFirst->GetXaxis()->SetTitle("sin(theta)");
 mygraphFirst->GetYaxis()->SetTitle("1/sqrt(V) (volt^(-1/2))");




mygraphSecond->Draw("A*");
TCanvas *c2 = new TCanvas();
TF1 *fsecond = new TF1("fsecond","[0]*x+[1]",0.036,0.052);
mygraphSecond->Fit(fsecond);
double slopeSecond = fsecond->GetParameter(0);
double slope_error_second = fsecond->GetParError(0);
double d_second = slopeSecond * h / sqrt(8 * q * m);
double d_error_second = slope_error_second * h / sqrt(8 *
q * m);
 mygraphSecond->GetXaxis()->SetTitle("sin(theta)");
 mygraphSecond->GetYaxis()->SetTitle("1/sqrt(V) (volt^(-1/2))");


gStyle->SetOptFit(1111);

cout<<"the lattice constant d for the first ring is d_first = "<<d_first<<"+-"<<d_error_first<<
" and for the second ring is d_second =  "<<d_second<<"+-"<<d_error_second<< "\n";











}