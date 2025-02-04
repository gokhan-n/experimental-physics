
{
// PART 0 
const int ndata0 = 9;
float erI0[ndata0] = {0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
float voltages0[ndata0] = {     
   0.0347,0.0438,0.050,0.0566,0.0632,0.0689,0.0778,0.0855,0.0954
    };

float currents0[ndata0] = {     
   0.109,0.137,0.156,0.176,0.195,0.211,0.237,0.258,0.285
    };

float erV0[ndata0];

for (int p = 0; p < ndata0; ++p) {
    if (voltages0[p] <= 0.4)
{
erV0[p] = (3 * voltages0[p] / 100) + 0.0004;
}
else if (voltages0[p] <= 4)
{
erV0[p] = (5 * voltages0[p] / 100) + 0.003;
}
else if (voltages0[p] <= 40)
{
erV0[p] = (5 * voltages0[p] / 100) + 0.03;
}
else if (voltages0[p] <= 400)
{
    erV0[p] = (5 * voltages0[p] / 100) + 0.3;
}
else
{
erV0[p] = (voltages0[p] / 100) + 4;
}

}

TGraphErrors *part0 = new TGraphErrors(ndata0 ,currents0,voltages0,erI0,erV0);
part0->Draw("A*");
part0->GetXaxis()->SetTitle("Current(A)");
part0->GetYaxis()->SetTitle("Voltage(V)");
part0->SetTitle("I vs V for resistence");


TF1 *p0Linefit = new TF1("p0Linefit","[0]*x+[1]");
p0Linefit->SetParameters(3.1416,2.7182); 
part0->Fit(p0Linefit);
TCanvas *c = new TCanvas();

gStyle->SetOptFit(1111);

double Rroom = p0Linefit->GetParameter(0);
double Rroom_error = p0Linefit->GetParError(0);

cout<<"resistance at room temperature 300K is  = "<<Rroom<<"+-"<<Rroom_error<< "\n";

// T VS R/R_300K

const int ntable = 34;
float Ttable[ntable] = {300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,
2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500};
float TtableEr[ntable] = {100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,
100,100,100,100,100,100,100,100,100,100,100};
float Rtable[ntable] = {1.0,1.43,1.87,2.34,2.85,3.36,3.88,4.41,4.95,5.48,6.03,6.58,7.14,7.71,8.28,8.86,9.44,10.03,10.63,11.24,11.84,12.46,
13.08,13.72,14.34,14.99,15.63,16.29,16.95,17.62,18.28,18.97,19.66};
float RtableEr[ntable] = {0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,
0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};

TGraphErrors *table = new TGraphErrors(ntable,Rtable,Ttable,RtableEr,TtableEr);
table->Draw("A*");
table->GetXaxis()->SetTitle("R/R_300K");
table->GetYaxis()->SetTitle("T(Kelvin)");
table->SetTitle("T vs R/R_300K given values");


TF1 *tablepoly = new TF1("tableLine","[0]*x*x + [1]*x + [2]");
tablepoly->SetParameters(1,1,1); 
table->Fit(tablepoly);
TCanvas *c2 = new TCanvas();

double tSecond = tablepoly->GetParameter(0);
double tFirst = tablepoly->GetParameter(1);
double tZero = tablepoly->GetParameter(2);

cout<<"T(R/R_300K) = "<<tSecond<<"(R/R_300K)^2"<<tFirst<<"R/R_300K"<<tZero<< "\n";







// PART 1

const int ndata1 = 14;
float amp = 1e4;
float vLamb1[ndata1] = {0.06,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.70,0.80};
float vLamb1err[ndata1] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};

float vLamb1log[ndata1];
float vLamb1logerr[ndata1];





float d1[ndata1] = {32.7,24.2,19.4,17.2,15.4,13.6,12.5,11.8,11.3,10.5,10.0,9.7,9.0,8.3};
float d1err[ndata1] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};

float d1log[ndata1];
float d1logerr[ndata1];

for (int i = 0; i < 14; ++i) {
    vLamb1[i] = vLamb1[i] / amp;
    vLamb1err[i] = vLamb1err[i] / amp;
    vLamb1log[i] = log(vLamb1[i]);
    vLamb1logerr[i] = vLamb1err[i] / vLamb1[i];
    d1log[i] = log(d1[i]);
    d1logerr[i] = d1err[i] / d1[i];

}



TGraphErrors *p1_exp = new TGraphErrors(ndata1,d1,vLamb1,d1err,vLamb1err);
p1_exp->Draw("A*");
p1_exp->GetXaxis()->SetTitle("distance(cm)");
p1_exp->GetYaxis()->SetTitle("Voltage(V)");
p1_exp->SetTitle("distance vs V- exp-fit ");


TF1 *exp = new TF1("exp","[1]*pow(x,[0])");
exp->SetParNames("Power", "Scale");
p1_exp->Fit(exp);

TCanvas *c3 = new TCanvas();



float n = exp->GetParameter(0);
float n_error = exp->GetParError(0);
cout << "Power = "<< n << "+- "<<n_error<<"\n";

TGraphErrors *p1log = new TGraphErrors(ndata1,d1log,vLamb1log,d1logerr,vLamb1logerr);
p1log->Draw("A*");
p1log->GetXaxis()->SetTitle("Log(distance(cm))");
p1log->GetYaxis()->SetTitle("Log(Voltage(V))");
p1log->SetTitle("distance vs V- log-fit ");


TF1 *p1Linefit = new TF1("p1Linefit","[0]*x+[1]");
p1log->Fit(p1Linefit);
TCanvas *c4 = new TCanvas();

// PART 2 

// data for cm


const int ndata2 = 10;
float p2Vsensor7cm[ndata2] = {1.0,0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60,0.50};
float p2Vlamb7cm[ndata2] = {3.932,3.721,3.661,3.417,3.246,3.114,2.953,2.808,2.648,2.245};
float p2Ilamb7cm[ndata2]= {1.689,1.644,1.631,1.579,1.542,1.513,1.477,1.444,1.407,1.311};

//data for 10 cm
float p2Vsensor10cm[9]= {1.0,0.95,0.90,0.85,0.80,0.75,0.70,0.65,0.60};
float p2Vlamb10cm[9]= {4.86,4.58,4.39,4.20,3.99,3.83,3.66,3.493,3.298};
float p2Ilamb10cm[9]= {1.875,1.818,1.782,1.742,1.700,1.665,1.630,1.594,1.553};

// floats to be filled

float p2Vsensor7cmerr[ndata2],p2Vlamb7cmerr[ndata2],p2Ilamb7cmerr[ndata2],
p2Rlamb7cm[ndata2],p2Rlamb7cmerr[ndata2],T7cm[ndata2],T7cmerr[ndata2]
;

float p2Vsensor10cmerr[9],p2Vlamb10cmerr[9],p2Ilamb10cmerr[9],
p2Rlamb10cm[9],p2Rlamb10cmerr[9],T10cm[9],T10cmerr[9]
;

// log floats

float T7cmlog[ndata2], T7cmlogerr[ndata2], T10cmlog[9],T10cmlogerr[9],
p2Vsensor7cmlog[ndata2],p2Vsensor7cmlogerr[ndata2], p2Vsensor10cmlog[9],p2Vsensor10cmlogerr[9]; 

// for 7 cm fgetting errors and R values

for (int j = 0; j < 10; ++j) {

    // sensor in terms of voltage
    p2Vsensor7cm[j] = p2Vsensor7cm[j] / amp;
    // errors of measurements
    p2Vsensor7cmerr[j] = 0.05 / amp;
    p2Vlamb7cmerr[j] = 0.001;
    p2Ilamb7cmerr[j] = 0.001;
    //resistance calculations
    p2Rlamb7cm[j] = p2Vlamb7cm[j] / p2Ilamb7cm[j];
    p2Rlamb7cmerr[j] = p2Rlamb7cm[j] * sqrt(pow(p2Vlamb7cmerr[j]/p2Vlamb7cm[j],2)+
    pow(p2Ilamb7cmerr[j]/p2Ilamb7cm[j],2));


    


    // getting T and its error
    
    T7cm[j] = tSecond * (p2Rlamb7cm[j]/Rroom)*(p2Rlamb7cm[j]/Rroom) + 
    tFirst * (p2Rlamb7cm[j]/Rroom) + tZero;

    double dTdR = tSecond * (2 * p2Rlamb7cm[j]) / (Rroom * Rroom) + tFirst / Rroom;
    double dTdRroom = -tSecond * (2 * p2Rlamb7cm[j] * p2Rlamb7cm[j]) / (Rroom * Rroom * Rroom) - tFirst * p2Rlamb7cm[j] / (Rroom * Rroom);

    T7cmerr[j] = sqrt((dTdR * p2Rlamb7cmerr[j]) * (dTdR * p2Rlamb7cmerr[j]) +
     (dTdRroom * Rroom_error) * (dTdRroom * Rroom_error));
     
    cout<<"for 7 cm T values are: "<<T7cm[j]<<"+-"<<T7cmerr[j]<<"\n";

    //filling logs

    p2Vsensor7cmlog[j] = log(p2Vsensor7cm[j]);
    p2Vsensor7cmlogerr[j] = p2Vsensor7cmerr[j] / p2Vsensor7cm[j];

    T7cmlog[j] = log(T7cm[j]);
    T7cmlogerr[j] = T7cmerr[j]/T7cm[j];
}

// for 10 cm getting errors and R values
for (int k = 0; k < 9; ++k) {

    // sensor in terms of voltage
    p2Vsensor10cm[k] = p2Vsensor10cm[k] / amp;
    // errors of measurements
    p2Vsensor10cmerr[k] = 0.05 / amp;
    p2Vlamb10cmerr[k] = 0.001;
    p2Ilamb10cmerr[k] = 0.001;
    //resistance calculations
    p2Rlamb10cm[k] = p2Vlamb10cm[k] / p2Ilamb10cm[k];
    p2Rlamb10cmerr[k] = p2Rlamb10cm[k] * sqrt(pow(p2Vlamb10cmerr[k]/p2Vlamb10cm[k],2)+
    pow(p2Ilamb10cmerr[k]/p2Ilamb10cm[k],2));
    


    // getting T and its error
    
    T10cm[k] = tSecond * (p2Rlamb10cm[k]/Rroom)*(p2Rlamb10cm[k]/Rroom) + 
    tFirst * (p2Rlamb10cm[k]/Rroom) + tZero;
    

    double dTdR_10cm = tSecond * (2 * p2Rlamb10cm[k]) / (Rroom * Rroom) + tFirst / Rroom;
    double dTdRroom_10cm = -tSecond * (2 * p2Rlamb10cm[k] * p2Rlamb10cm[k]) / (Rroom * Rroom * Rroom) - tFirst * p2Rlamb10cm[k] / (Rroom * Rroom);

    T10cmerr[k] = sqrt((dTdR_10cm * p2Rlamb10cmerr[k]) * (dTdR_10cm * p2Rlamb10cmerr[k]) +
     (dTdRroom_10cm * Rroom_error) * (dTdRroom_10cm * Rroom_error));
     

    // filling logs 

    p2Vsensor10cmlog[k] = log(p2Vsensor10cm[k]);
    p2Vsensor10cmlogerr[k] = p2Vsensor10cmerr[k] / p2Vsensor10cm[k];

    T10cmlog[k] = log(T10cm[k]);
    T10cmlogerr[k] = T10cmerr[k]/T10cm[k];

    cout<<"for 7 cm T values are: "<<T10cm[k]<<"+-"<<T10cmerr[k]<<"\n";
    
}

    
// graphs for measurements of 7 cm distance

 // exponantial fit
    
    TGraphErrors *p2exp_7cm = new TGraphErrors(ndata2,T7cm,p2Vsensor7cm,T7cmerr,p2Vsensor7cmerr);
    p2exp_7cm->Draw("A*");
    p2exp_7cm->GetXaxis()->SetTitle("T(Kelvin)");
    p2exp_7cm->GetYaxis()->SetTitle("Voltage(V)");
    p2exp_7cm->SetTitle("T vs Vsensor for 7 cm distance");

    
    TF1 *exp_7cm = new TF1("exp7","[1]*TMath::Power(x,[0])");
    exp_7cm->SetParameter(1, 1.70743068e-17);
    exp_7cm->SetParameter(0, 4);
    p2exp_7cm->Fit(exp_7cm);
    exp_7cm->SetParNames("Power", "Scale");


    TCanvas *c5 = new TCanvas();

// log fit

    TGraphErrors *p2log_7cm = new TGraphErrors(ndata2,T7cmlog,p2Vsensor7cmlog,T7cmlogerr,p2Vsensor7cmlogerr);
    p2log_7cm->Draw("A*");
    p2log_7cm->GetXaxis()->SetTitle("log(T(Kelvin))");
    p2log_7cm->GetYaxis()->SetTitle("log(Voltage(V))");
    p2log_7cm->SetTitle("T vs Vsensor for 7 cm distance -logfit");

    
    TF1 *log_7cm = new TF1("log7","[0]*x+[1]");
    p2log_7cm->Fit(log_7cm);
    log_7cm->SetParNames("slope", "intercept");


    TCanvas *c6= new TCanvas();


    // graphs for measurements of 10 cm distance

 // exponantial fit

  const int n10cm = 9;
    
    TGraphErrors *p2exp_10cm = new TGraphErrors(n10cm,T10cm,p2Vsensor10cm,T10cmerr,p2Vsensor10cmerr);
    p2exp_10cm->Draw("A*");
    p2exp_10cm->GetXaxis()->SetTitle("T(Kelvin)");
    p2exp_10cm->GetYaxis()->SetTitle("Voltage(V)");
    p2exp_10cm->SetTitle("T vs Vsensor for 10 cm distance");
    
    TF1 *exp_10cm = new TF1("exp10","[1]*pow(x,[0])");
    exp_10cm->SetParameter(1, 1.70743068e-17);
    exp_10cm->SetParameter(0, 4);
    p2exp_10cm->Fit(exp_10cm);
    exp_10cm->SetParNames("Power", "Scale");

    TCanvas *c7 = new TCanvas();
    

    // log fit

    TGraphErrors *p2log_10cm = new TGraphErrors(n10cm,T10cmlog,p2Vsensor10cmlog,T10cmlogerr,p2Vsensor10cmlogerr);
    p2log_10cm->Draw("A*");
    p2log_10cm->GetXaxis()->SetTitle("log(T(Kelvin))");
    p2log_10cm->GetYaxis()->SetTitle("log(Voltage(V))");
    p2log_10cm->SetTitle("T vs Vsensor  for 10 cm distance -logfit");

    
    TF1 *log_10cm = new TF1("log10","[0]*x+[1]");
    p2log_10cm->Fit(log_10cm);
    log_10cm->SetParNames("slope", "intercept");


    TCanvas *c8= new TCanvas();



// PART 3

const int ndata3 =20;
float p3T[ndata3] = {250.6,252.1,255.3,258.6,261.1,263.0,269.4,272.9,278.3,282.6,286.9,291.3,295.6,299.8,
304.1,308.2,312.0,316.1,320.0,324.0};
float p3Vsensor[ndata3] = {4.8,4.9,5.0,5.1,5.2,5.3,5.5,5.7,5.9,6.1,6.3,
6.5,6.7,6.9,7.1,7.3,7.5,7.7,7.9,8.1};

float p3Terr[ndata3],p3Vsensorerr[ndata3],p3Tlog[ndata3],p3Vlog[ndata3],p3Tlogerr[ndata3],p3Vlogerr[ndata3];

for (int m = 0; m < ndata3; ++m) {

p3T[m] = p3T[m] + 273.1;
p3Terr[m] = 0.1;
p3Vsensor[m] = p3Vsensor[m] / 1000;
p3Vsensorerr[m] = 0.1 /1000;
p3Tlog[m] = log(p3T[m]);
p3Vlog[m] = log(p3Vsensor[m]);


p3Tlogerr[m] = p3Terr[m]/p3T[m];
p3Vlogerr[m] = p3Vsensorerr[m]/p3Vsensor[m];

}

// part 3 exponantial
    TGraphErrors *p3exp = new TGraphErrors(ndata3,p3T,p3Vsensor,p3Terr,p3Vsensorerr);
    p3exp->Draw("A*");
    p3exp->GetXaxis()->SetTitle("T(Kelvin)");
    p3exp->GetYaxis()->SetTitle("Voltage(V)");
    p3exp->SetTitle("Low Temperature Exp-fit");
    
    TF1 *lastexp = new TF1("lastexp","[1]*pow(x,[0])");
    lastexp->SetParameter(1, 1.70743068e-17);
    lastexp->SetParameter(0, 4);
    p3exp->Fit(lastexp);
    lastexp->SetParNames("Power", "Scale");

    TCanvas *c9 = new TCanvas();

// part 3 log

TGraphErrors *p3log = new TGraphErrors(ndata3,p3Tlog,p3Vlog,p3Tlogerr,p3Vlogerr);
    p3log->Draw("A*");
    p3log->GetXaxis()->SetTitle("log(T(Kelvin))");
    p3log->GetYaxis()->SetTitle("log(Voltage(V))");
    p3log->SetTitle("Low Temperature Log-fit");
    
    TF1 *lastlog = new TF1("lastlog","[0]*x+[1]");
    p3log->Fit(lastlog);
    lastlog->SetParNames("slope", "intercept");

    TCanvas *c10 = new TCanvas();










}


















































































































