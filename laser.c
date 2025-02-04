{
    //SINGLE&DOUBLE SLIT  VALUES 

    float diffLaser = 650*1e-9;
    float diffLaser_error = 0.1*1e-7;
    float L_double = 1.511;
    float L_single = 1.497;
    float L_hair = 1.510;
    float L_double_error = 0.001;
    float L_single_error = 0.001;
    float L_hair_error = 0.001;
    float d_25 = 0.25e-3;
    float d_25_error = 0.01e-3;
    float d_50 = 0.50e-3;
    float d_50_error = 0.01e-3;
    float a = 0.16e-3;
    float a_error = 0.01e-3;

    // MICHELSON  VALUES

    float michLaser= 632.8*1e-9;
    float michLaser_error = 0.1*1e-7;
    float degree_error = 0.1;
    float d_mich = 0.55*1e-2;
    float d_mich_error = 0.001*1e-2;

    // PFUND VALUES 

    float r = 1.0*1e-2;
    float r_error = 0.1*1e-2;



    //DATA
    float d25Left[4]={0.4,0.8,1.1,1.5};
    float d25Right[4]={0.4,0.8,1.2,1.5};

    float d50Left[4]={0.2,0.4,0.5,0.7};
    float d50Right[4]={0.2,0.4,0.6,0.8};


    float sLeft[4]={0.7,1.4,1.9,2.5};
    float sRight[4]={0.6,1.2,1.8,2.3};

    float sUpper[4]={1.3,2.8,4.2,5.5};
    float sLower[4]={1.0,2.8,4.2,5.6};

    float mValues[10]={10,20,30,40,50,60,70,80,90,100};
    float thetas[10]={3.3,4.7,5.8,6.8,7.5,8.2,8.9,9.4,9.9,10.4};

    float pi = TMath::Pi();
   for (int i = 0; i < 10; ++i) {
       thetas[i] = thetas[i]*pi/180;
    }
    


    //CALCULATIONS FOR SINGLE&DOUBLE SLITS

    float d25[4],d50[4],single[4],thickness[4],nvalues[10],nvalues_errors[10];

    for (int i = 0; i < 4; ++i) {
        d25[i] =1e-2*(d25Left[i]+d25Right[i])/2;
        d50[i] = 1e-2*(d50Left[i]+d50Right[i])/2;
        single[i] = 1e-2*(sLeft[i]+sRight[i])/2;
        thickness[i] = 1e-2*(sUpper[i]+sLower[i])/2;
    }

    float d25Lambdas[4],d50Lambdas[4],sWidth[4],hairWidth[4];
    float d25Lambdas_errors[4],d50Lambdas_errors[4],sWidth_errors[4],hairWidth_errors[4];

    // DOUBLE SLIT D25
    for (int i = 0; i < 4; ++i) {
    int n = i + 1;
    float Delta_y = d25[i];
    float sqrt_term = sqrt(pow(L_double, 2) + pow(Delta_y, 2));

    // Calculate d25Lambdas
    d25Lambdas[i] = 1e9 * d_25 * Delta_y / (n * sqrt_term);

    // Partial derivatives
    float partial_d = Delta_y / (n * sqrt_term);
    float partial_Delta_y = (d_25 / (n * sqrt_term)) - (d_25 * Delta_y / (n * pow(sqrt_term, 3))) * Delta_y;
    float partial_L = -(d_25 * Delta_y * L_double) / (n * pow(sqrt_term, 3));

    // Error contributions
    float error_d = partial_d * d_25_error;
    float error_Delta_y = partial_Delta_y * (1e-2 * 0.01); // Assuming 0.01 cm error in Delta_y
    float error_L = partial_L * L_double_error;

    // Total error
    d25Lambdas_errors[i] = 1e9 * sqrt(pow(error_d, 2) + pow(error_Delta_y, 2) + pow(error_L, 2));

    // Print results
    cout << "Lambda for d25[" << i << "] = " << d25Lambdas[i] << " ± " << d25Lambdas_errors[i] << " nm" << endl;
}

    //DOUBLE SLIT D50
    for (int i = 0; i < 4; ++i) {
    int n = i + 1;
    float Delta_y = d50[i];
    float sqrt_term = sqrt(pow(L_double, 2) + pow(Delta_y, 2));

    // Calculate d50Lambdas
    d50Lambdas[i] = 1e9 * d_50 * Delta_y / (n * sqrt_term);

    // Partial derivatives
    float partial_d = Delta_y / (n * sqrt_term);
    float partial_Delta_y = (d_50 / (n * sqrt_term)) - (d_50 * Delta_y / (n * pow(sqrt_term, 3))) * Delta_y;
    float partial_L = -(d_50 * Delta_y * L_double) / (n * pow(sqrt_term, 3));

    // Error contributions
    float error_d = partial_d * d_50_error;
    float error_Delta_y = partial_Delta_y * (1e-2 * 0.01); // Assuming 0.01 cm error in Delta_y
    float error_L = partial_L * L_double_error;

    // Total error
    d50Lambdas_errors[i] = 1e9 * sqrt(pow(error_d, 2) + pow(error_Delta_y, 2) + pow(error_L, 2));

    // Print results
    cout << "Lambda for d50[" << i << "] = " << d50Lambdas[i] << " ± " << d50Lambdas_errors[i] << " nm" << endl;
}


 
    // SINGLE SLIT
   for (int i = 0; i < 4; ++i) {
    int n = i + 1;
    float Delta_y = single[i];
    float sqrt_term = sqrt(pow(L_single, 2) + pow(Delta_y, 2));

    // Calculate single slit thickness (sWidth)
    sWidth[i] = 1e3 * n * diffLaser * (n * sqrt_term) / Delta_y;

    // Partial derivatives
    float partial_Delta_y = -1e3 * n * diffLaser * (n * sqrt_term / pow(Delta_y, 2)) +
                            1e3 * n * diffLaser * (n * Delta_y / (sqrt_term * Delta_y));
    float partial_L = 1e3 * n * diffLaser * (n * L_single / (sqrt_term * Delta_y));
    float partial_diffLaser = 1e3 * n * (n * sqrt_term) / Delta_y;

    // Error contributions
    float error_Delta_y = partial_Delta_y * (1e-2 * 0.01); // Assuming 0.01 cm error in Delta_y
    float error_L = partial_L * L_single_error;
    float error_diffLaser = partial_diffLaser * diffLaser_error;

    // Total error
    sWidth_errors[i] = sqrt(pow(error_Delta_y, 2) + pow(error_L, 2) + pow(error_diffLaser, 2));

    // Print results
    cout << "Single slit thickness for order " << n << " = " << sWidth[i] << " ± " << sWidth_errors[i] << " mm" << endl;
}
    // SINGLE SLIT THICKNESS
    for (int i = 0; i < 4; ++i) {
    int n = i + 1;
    float Delta_y = thickness[i];
    float sqrt_term = sqrt(pow(L_hair, 2) + pow(Delta_y, 2));

    // Calculate hair thickness (hairWidth)
    hairWidth[i] = 1e3 * n * diffLaser * (n * sqrt_term) / Delta_y;

    // Partial derivatives
    float partial_Delta_y = -1e3 * n * diffLaser * (n * sqrt_term / pow(Delta_y, 2)) +
                            1e3 * n * diffLaser * (n * Delta_y / (sqrt_term * Delta_y));
    float partial_L = 1e3 * n * diffLaser * (n * L_hair / (sqrt_term * Delta_y));
    float partial_diffLaser = 1e3 * n * (n * sqrt_term) / Delta_y;

    // Error contributions
    float error_Delta_y = partial_Delta_y * (1e-2 * 0.01); // Assuming 0.01 cm error in Delta_y
    float error_L = partial_L * L_hair_error;
    float error_diffLaser = partial_diffLaser * diffLaser_error;

    // Total error
    hairWidth_errors[i] = sqrt(pow(error_Delta_y, 2) + pow(error_L, 2) + pow(error_diffLaser, 2));

    // Print results
    cout << "Hair thickness for order " << n << " = " << hairWidth[i] << " ± " << hairWidth_errors[i] << " mm" << endl;
}

   

    // CALCULATIONS FOR MICHELSON
    for (int i = 0; i < 10; ++i) {
    float m = mValues[i];
    float cos_theta = cos(thetas[i]);
    float one_minus_cos = 1 - cos_theta;

    // Calculate index of refraction
    nvalues[i] = (d_mich - (m * michLaser / 2)) * one_minus_cos /
                 (d_mich * one_minus_cos - (m * michLaser / 2));

    // Partial derivatives
    float partial_d_mich = one_minus_cos / (d_mich * one_minus_cos - (m * michLaser / 2)) -
                           ((d_mich - (m * michLaser / 2)) * one_minus_cos) /
                           pow(d_mich * one_minus_cos - (m * michLaser / 2), 2);

    float partial_michLaser = -m / 2 * one_minus_cos /
                              (d_mich * one_minus_cos - (m * michLaser / 2)) +
                              ((d_mich - (m * michLaser / 2)) * m / 2 * one_minus_cos) /
                              pow(d_mich * one_minus_cos - (m * michLaser / 2), 2);

    float partial_theta = (d_mich - (m * michLaser / 2)) * sin(thetas[i]) /
                          (d_mich * one_minus_cos - (m * michLaser / 2)) -
                          ((d_mich - (m * michLaser / 2)) * one_minus_cos * d_mich * sin(thetas[i])) /
                          pow(d_mich * one_minus_cos - (m * michLaser / 2), 2);

    // Error contributions
    float error_d_mich = partial_d_mich * d_mich_error;
    float error_michLaser = partial_michLaser * michLaser_error;
    float error_theta = partial_theta * degree_error * TMath::Pi() / 180; // Convert degree error to radians

    // Total error
    nvalues_errors[i] = sqrt(pow(error_d_mich, 2) + pow(error_michLaser, 2) + pow(error_theta, 2));

    // Print results
    cout << "Index of refraction for m = " << m << ": " << nvalues[i]
         << " ± " << nvalues_errors[i] << endl;
}

  

    // CALCULATIONS FOR PFUND

    // Calculate Pfund index of refraction
    float nPfund = sqrt(1 + pow(2 * d_mich / r, 2));

// Partial derivatives
float partial_d_mich = (4 * d_mich / pow(r, 2)) / sqrt(1 + pow(2 * d_mich / r, 2));
float partial_r = (-4 * pow(d_mich, 2) / pow(r, 3)) / sqrt(1 + pow(2 * d_mich / r, 2));

// Error contributions
float error_d_mich = partial_d_mich * d_mich_error;
float error_r = partial_r * r_error;

// Total error
float nPfund_error = sqrt(pow(error_d_mich, 2) + pow(error_r, 2));

// Print result
cout << "Index of refraction for Pfund = " << nPfund << " ± " << nPfund_error << endl;


    



}