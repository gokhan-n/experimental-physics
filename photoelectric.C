{
  std::vector<std::string> voltages = {
        "yellow-V.txt",
        "green-V.txt",
        "turquoise-V.txt",
        "blue-V.txt",
        "violet-V.txt"


     };

std::vector<std::string> currents = {
        "yellow-I.txt",
        "green-I.txt",
        "turquoise-I.txt",
        "blue-I.txt",
        "violet-I.txt"


};std::vector<std::string> names = {
        "yellow",
        "green",
        "turquoise",
        "blue",
        "violet"


};

std::vector<std::vector<double>> fitRangesFirst = {
      

    {0,0.200},{0,0.370},{0,0.330},{0,0.260},{0,0.727}
          
    };

    /*std::vector<std::vector<double>> fitRangesFirst = {
      

    {0.025,0.220},{0.205,0.370},{0.145,0.330},{0.065,0.230},{0.574,0.850}
          
    };
    */

std::vector<std::vector<double>> fitRangesSecond = {
      

    {2.500,3.000},{2.500,3.000},{2.150,3.000},{2.100,3.000},{2.010,2.980}
          
    };
std::vector<std::vector<double>> parameterRanges = {
      

    {0.4,-1.0},{-12.39,5.66},{-12,4.7},{-12.21,5.264},{-38.27,32.18}
          
    };

 
float del_I = 0.01;
float del_V = 0.1*pow(10,-3);
float freq [5]={519e12,549e12,608e12,688e12,711e12};
float stopping[5];
float error_sy[5];



for(int i =0;i<voltages.size();i++){
    
    
    std::ifstream fileVoltage(voltages[i].c_str());
    int ndata;

    std::string line;
    int linecount = 0; // Initialize a counter
    while (std::getline(fileVoltage, line)) {  // Loop through each line in the file
        linecount++;  // Incrementing line count for each line read

    }

    ndata = linecount;
    float x[ndata], y[ndata],sx[ndata],sy[ndata];
    std::vector<float> retarded;
    std::vector<float> amplified;
   
    std::ifstream infile1(voltages[i]);
    float value1;

    while (infile1 >> value1) {
        float InsertV = value1*pow(10,-3);
        retarded.push_back(value1);
    }

    std::ifstream infile2(currents[i]);
    float value2;

    while (infile2 >> value2) {
        float InsertI = value2*pow(10,-12);
        amplified.push_back(InsertI);
    }

    for(int k =0;k<ndata;k++){
        x[k] = retarded[k]/pow(10,3);
        y[k]= amplified[k]*pow(10,12);
        sx[k] = del_V;
        if(i<3){
            sy[k] = 0.001;
        }
        else {
            sy[k] = del_I;

        }
    }
       
       


     /*for(int p = 0;p<ndata;p++){
         cout<< p<<"th voltage element is"<<x[p]<<"\n";

     }
     for(int r = 0;r<ndata;r++){
         cout<< r<<"th current element is"<<y[r]<<"\n";

     }

     for(int t = 0;t<ndata;t++){
         cout<< t<<"th elements are"<<x[t]<<" V "<<y[t]<<" I "<<"\n";

     }
     */

    TGraphErrors *mygraph = new TGraphErrors(ndata,x,y,sx,sy);
    TCanvas *c = new TCanvas();


    mygraph->Draw("A*");
    mygraph->SetTitle(names[i].c_str());
    mygraph->GetXaxis()->SetTitle("Voltage[V]");
    mygraph->GetYaxis()->SetTitle("Current[10^-12 A]");
   
    TF1 *line1 = new TF1("line","[0]*x+[1]",fitRangesFirst[i][0],fitRangesFirst[i][1]);
    line1->SetParameter(parameterRanges[i][0],parameterRanges[i][1]);
    mygraph->Fit(line1,"R");

    TF1 *line2 = new TF1("line","[0]*x+[1]",fitRangesSecond[i][0],fitRangesSecond[i][1]);
    mygraph->Fit(line2,"R+");

    float m1 = line1->GetParameter(1);
    float n1 = line1->GetParameter(0);
    float m2 = line2->GetParameter(1);
    float n2 = line2->GetParameter(0);

    float m1_error = line1->GetParError(1);
    float n1_error = line1->GetParError(0);
    float m2_error = line2->GetParError(1);
    float n2_error = line2->GetParError(0);
    
    cout<<"for m1,n1,m2,n2: "<<m1<<","<<n1<<","<<m2<<","<<n2<<" errors are "
    <<m1_error<<","<<n1_error<<","<<m2_error<<","<<n2_error<<endl;

    cout<<"for"<< currents[i].c_str()<<endl;
    cout<<"for first line "<<"\n"; 
    cout <<"intercept is "<< n1 <<"+-"<<n1_error<<"slope is "<<m1<<"+-"<<m1_error<<"\n";
    cout<<"for second line "<<"\n"; 
    cout <<"intercept is "<<n2 <<"+-"<<n2_error<<"slope is "<< m2<<"+-"<<m2_error<<"\n";

    float Vs = (m2-m1)/(n1-n2);
    cout<<"Vs potential is "<<Vs<<endl;



    gStyle->SetOptFit(1111);

    

    stopping[i]= Vs;

    float up= m2-m1;
    float down= n1-n2;
    float up_error = sqrt(m1_error*m1_error+m2_error*m2_error);
    float down_error = sqrt(n1_error*n1_error+n2_error*n2_error);


    float error = Vs * sqrt((up_error*up_error)/(up*up)+(down_error*down_error)/(down*down));

    error_sy[i]=error;



    


    cout<< "final Vs value for "<<names[i].c_str()<<"is " << Vs <<"+-"<<error<< endl;

    
}
    TGraphErrors *finalgraph = new TGraphErrors(5,freq,stopping,0,error_sy);
    TCanvas *c6 = new TCanvas();
    finalgraph->Draw("A*");
    TF1 *finalFit = new TF1("finalFit","[0]*x+[1]");
    finalgraph->Fit(finalFit);
    finalgraph->SetTitle("Stopping voltage vs frequency");
    finalgraph->GetXaxis()->SetTitle("Frequency(Hz)");
    finalgraph->GetYaxis()->SetTitle("Voltage[V]");
    float m = finalFit->GetParameter(1);
    float n = finalFit->GetParameter(0);
    float m_error = finalFit->GetParError(1);
    float n_error = finalFit->GetParError(0);

    cout<<"the final value is"<<n<<"+-"<<n_error<<"with the y-intercept "<<m<<"+-"<<m_error<<endl;


}




