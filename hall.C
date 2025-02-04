{
double calculateVoltUncertainty(double V)
{
V = abs(V);
double error;
if (V <= 0.4)
{
error = (3 * V / 100);
}
else if (V <= 4)
{
error = (5 * V / 100) + 0.003;
}
else if (V <= 40)
{
error = (5 * V / 100) + 0.03;
}
else if (V <= 400)
{
error = (5 * V / 100) + 0.3;
}
else
{
error = (V / 100) + 4;
}
return error;
};

double t = 1e-3;
double t_error = 0.2e-3;
double L = 16e-3;
double L_error = 4e-3;
double w = 10e-3;
double w_error = 0.2e-3;
double A = w*t;
double A_error = sqrt(pow(w_error*t,2)+pow(w
*t_error,2));

// PART 1

const int ndata =10;
float I1[ndata] = {0.020,0.026,0.037,0.045,0.053,0.063,0.060,0.074,0.083,0.094};
float V1[ndata] = {0.907,1.154,1.650,2.031,2.369,2.827,3.126,3.366,3.765,4.280};
float V1er[ndata],I1er[ndata];

for (int i = 0; i < ndata; ++i) {
    I1er[i] = 0.001;
    V1er[i] = calculateVoltUncertainty(V1[i]); 
}

TGraphErrors *p1 = new TGraphErrors(ndata,I1,V1,I1er,V1er);
p1->Draw("A*");
p1->GetXaxis()->SetTitle("Current(A)");
p1->GetYaxis()->SetTitle("Voltage(V)");
p1->SetTitle("Finding Resistance-Linefit ");


TF1 *linefit = new TF1("linefit","[0]*x+[1]");
linefit->SetParNames("Slope", "Intercept");
p1->Fit(linefit);
gStyle->SetOptFit(1);

TCanvas *c = new TCanvas();

double R = linefit->GetParameter(0);
double Rerror = linefit->GetParError(0);

// for mobility calculations 
double sigma = L/(A*R);
double sigma_error = sqrt(pow(L_error/(A*R)
,2)+pow(A_error*sigma/A,2)+pow(Rerror*
sigma/R,2));
//

// PART 2   

// B = 0
float I2[ndata] = {0.012,0.019,0.026,0.032,0.047,0.056,0.063,0.069,0.075,0.083};
float V2[ndata] = {-0.00081,-0.00133,-0.00182,-0.00218,-0.00327,-0.00392,-0.00436,-0.00486,-0.00531,-0.00585};
float I2er[ndata],V2er[ndata];

for (int i = 0; i < ndata; ++i) {
    I2er[i] = 0.001;
    V2er[i] = calculateVoltUncertainty(V2[i]); 
}

TGraphErrors * nonmagnetic = new TGraphErrors(ndata,I2,V2,I2er,V2er);
nonmagnetic->Draw("A*");
nonmagnetic->GetXaxis()->SetTitle("Current(A)");
nonmagnetic->GetYaxis()->SetTitle("Voltage(V)");
nonmagnetic->SetTitle("Reference Measurements");

linefit->SetParameters(-1.34,0); 

nonmagnetic->Fit(linefit);
double offset = linefit->GetParameter(0);
double offset_constant = linefit->GetParameter(1);

cout<<"offset is "<< offset<<"+"<<offset_constant<<"\n";



TCanvas *c2 = new TCanvas();


// MAGNETIC 

double B = 139.3*1e-3;
double B_error = B *5/100;

// finding V/I and its error

double p2_second = t/B;
double p2_seconderr = (t/B)*sqrt(pow((t_error/t),2) + pow((B_error/B),2));


float I2mag[ndata] = {0.010,0.017,0.026,0.033,0.039,0.046,0.054,0.061,0.070,0.076};
float V2mag[ndata] = {-0.00183,-0.00296,-0.00453,-0.00589,-0.00690,-0.00802,-0.00956,-0.01073,-0.01233,-0.01382};
float I2mager[ndata],V2mager[ndata];

for (int i = 0; i < ndata; ++i) {
    I2mager[i] = 0.001;
    V2mager[i] = calculateVoltUncertainty(V2[i]); 
}

TGraphErrors * magnetic = new TGraphErrors(ndata,I2mag,V2mag,I2mager,V2mager);
magnetic->Draw("A*");
magnetic->GetXaxis()->SetTitle("Current(A)");
magnetic->GetYaxis()->SetTitle("Voltage(V)");
magnetic->SetTitle("I vs V_hall Measurements at B= 139.3 mT");

linefit->SetParameters(-1.34,0); 

magnetic->Fit(linefit);

TCanvas *c3 = new TCanvas();

//getting slope 
double p2_slope = linefit->GetParameter(0);
double p2_slope_error = linefit->GetParError(0);

// finding RH and its error 

double RH1 = p2_slope*p2_second;
double RH1_error = sqrt(pow(p2_slope,2)*pow(p2_seconderr,2)
+ pow(p2_second,2)*pow(p2_slope_error,2));

cout <<"RH1 is "<<RH1<<"+-"<<RH1_error<<"\n"; 

// mobility

double mu1 = RH1*sigma;
double mu1_error = sqrt(pow(RH1,2)*pow(sigma_error,2)
+ pow(sigma,2)*pow(RH1_error,2));

cout <<"mu1 is "<<mu1<<"+-"<<mu1_error<<"\n"; 



// PART 3 

// I = 53mA

double I3p = 53*1e-3;
double I3p_error = 1e-3;

// finding t/I and its error

double p3P_second = t/I3p;
double p3P_seconderr = (t/I3p)*sqrt(pow((t_error/t),2) + pow((I3p_error/I3p),2));


float Bp[ndata] = {0.0144,0.0264,0.0349,0.0473,0.0587,0.0707,0.0792,0.0932,0.1041,0.1268};
float V3p[ndata] = {-0.0431,-0.0479,-0.0514,-0.0563,-0.0610,-0.0658,-0.0692,-0.0749,-0.0793,-0.0884};
float Bper[ndata],V3per[ndata];

for (int i = 0; i < ndata; ++i) {
    Bper[i] = Bp[i]*5/100;
    V3per[i] = calculateVoltUncertainty(V3p[i]); 
}

TGraphErrors * positive = new TGraphErrors(ndata,Bp,V3p,Bper,V3per);
positive->Draw("A*");
positive->GetXaxis()->SetTitle("Magnetic Field(T)");
positive->GetYaxis()->SetTitle("Voltage(V)");
positive->SetTitle("B vs V_hall Measurements at I= 53mA");

linefit->SetParameters(-1.34,0); 

positive->Fit(linefit);

TCanvas *c4 = new TCanvas();

//getting slope 
double p3P_slope = linefit->GetParameter(0);
double p3P_slope_error = linefit->GetParError(0);

// finding RH and its error 

double RH2 = p3P_slope*p3P_second;
double RH2_error = sqrt(pow(p3P_slope,2)*pow(p3P_seconderr,2)
+ pow(p3P_second,2)*pow(p3P_slope_error,2));

cout <<"RH2 is "<<RH2<<"+-"<<RH2_error<<"\n"; 

// mobility 

double mu2 = RH2*sigma;
double mu2_error = sqrt(pow(RH2,2)*pow(sigma_error,2)
+ pow(sigma,2)*pow(RH2_error,2));

cout <<"mu2 is "<<mu2<<"+-"<<mu2_error<<"\n"; 



// I = -53mA
double I3n = -53*1e-3;
double I3n_error = 1e-3;

// finding t/I and its error

double p3N_second = t/I3n;
double p3N_seconderr = (t/I3n)*sqrt(pow((t_error/t),2) + pow((I3n_error/I3n),2));


float Bn[ndata] = {0.0187,0.0269,0.0332,0.0413,0.0499,0.0612,0.0712,0.0840,0.0975,0.1124};
float V3n[ndata] = {0.0448,0.0482,0.0507,0.0540,0.0574,0.0621,0.0661,0.0712,0.0767,0.0827};
float Bner[ndata],V3ner[ndata];

for (int i = 0; i < ndata; ++i) {
    Bner[i] = Bn[i]*5/100;
    V3ner[i] = calculateVoltUncertainty(V3n[i]); 
}

TGraphErrors * negative = new TGraphErrors(ndata,Bn,V3n,Bner,V3ner);
negative->Draw("A*");
negative->GetXaxis()->SetTitle("Magnetic Field(T)");
negative->GetYaxis()->SetTitle("Voltage(V)");
negative->SetTitle("B vs V_hall Measurements at I= -53mA");

linefit->SetParameters(1.34,0); 

negative->Fit(linefit);

TCanvas *c5 = new TCanvas();

//getting slope 
double p3N_slope = linefit->GetParameter(0);
double p3N_slope_error = linefit->GetParError(0);

// finding RH and its error 

double RH3 = p3N_slope*p3N_second;
double RH3_error = sqrt(pow(p3N_slope,2)*pow(p3N_seconderr,2)
+ pow(p3N_second,2)*pow(p3N_slope_error,2));

cout <<"RH3 is "<<RH3<<"+-"<<RH3_error<<"\n"; 

// mobility 

double mu3 = RH3*sigma;
double mu3_error = sqrt(pow(RH3,2)*pow(sigma_error,2)
+ pow(sigma,2)*pow(RH3_error,2));

cout <<"mu3 is "<<mu3<<"+-"<<mu3_error<<"\n"; 


// WITH OFFSET


float Voff[ndata],Voffer[ndata];

for (int i = 0; i < ndata; ++i) {
    Voff[i] = V2mag[i]- (offset*I2mag[i]+offset_constant);
    Voffer[i] = calculateVoltUncertainty(V2mag[i]- (offset*I2mag[i]+offset_constant));
 
}


TGraphErrors * withoff = new TGraphErrors(ndata,I2,Voff,I2er,Voffer);
withoff->Draw("A*");
withoff->GetXaxis()->SetTitle("Current(A)");
withoff->GetYaxis()->SetTitle("Voltage(V)");
withoff->SetTitle("The Measurements with Offset");

linefit->SetParameters(-1.34,0); 

withoff->Fit(linefit);

TCanvas *c6 = new TCanvas();

//getting slope 
double offset_slope = linefit->GetParameter(0);
double offset_slope_error = linefit->GetParError(0);

// finding RH and its error 

double RH4 = offset_slope*p2_second;
double RH4_error = sqrt(pow(offset_slope,2)*pow(p2_seconderr,2)
+ pow(p2_second,2)*pow(offset_slope_error,2));

cout <<"RH4 is "<<RH4<<"+-"<<RH4_error<<"\n"; 

// mobility 

double mu4 = RH4*sigma;
double mu4_error = sqrt(pow(RH4,2)*pow(sigma_error,2)
+ pow(sigma,2)*pow(RH4_error,2));

cout <<"mu4 is "<<mu4<<"+-"<<mu4_error<<"\n"; 

double RH_weighted = (RH1 / pow(RH1_error, 2) + RH2 / pow(RH2_error, 2) + RH3 / pow(RH3_error, 2) + RH4 / pow(RH4_error, 2)) /
                     (1 / pow(RH1_error, 2) + 1 / pow(RH2_error, 2) + 1 / pow(RH3_error, 2) + 1 / pow(RH4_error, 2));

double RH_weighted_error = sqrt(1/(1/pow(RH1_error,2)+1/pow(RH2_error,2)+1/pow(RH3_error,2)+1/pow(RH4_error,2)));

double mu_weighted = (mu1 / pow(mu1_error, 2) + mu2 / pow(mu2_error, 2) + mu3 / pow(mu3_error, 2) + mu4 / pow(mu4_error, 2)) /
                     (1 / pow(mu1_error, 2) + 1 / pow(mu2_error, 2) + 1 / pow(mu3_error, 2) + 1 / pow(mu4_error, 2));

double mu_weighted_error = sqrt(1/(1/pow(mu1_error,2)+1/pow(mu2_error,2)+1/pow(mu3_error,2)+1/pow(mu4_error,2)));



cout <<"RH weighted is "<<RH_weighted<<"+-"<<RH_weighted_error<<"\n"; 
cout <<"mu weighted is "<<mu_weighted<<"+-"<<mu_weighted_error<<"\n"; 


































}