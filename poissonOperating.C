 {
float operatingVoltages[4] = {380,400,420,440};
  float firstCounts[4] = {3922,3941,3982,3947};

 

  TGraphErrors *operating = new TGraphErrors(4,operatingVoltages,firstCounts,0,0);
  TCanvas *c1 = new TCanvas();

  operating->Draw("A*");
  operating->SetTitle("Voltage vs Counts;Volt(V);Counts");
  operating->GetXaxis()->SetLimits(300, 450);
  operating->SetMinimum(0);   
  operating->SetMaximum(4500); 

  TF1 *fnew = new TF1("fnew","[0]*x+[1]",360,500);
  fnew->SetParameters(0.01,1300);
  operating->Fit(fnew,"R");

}
