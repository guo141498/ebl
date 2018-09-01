#include<iostream>
#include<vector>
#include<nlopt.hpp>
#include<fstream>
#include<stdlib.h>
#include<map>
#include<cmath>
#include<random>
#include<complex.h>
#include "matrix.h"
#include <cstdlib>
#include<algorithm>
#include<time.h>
//version: tur
//ciber 
//parabola 
using namespace std;
double pi = 3.1415926;

//----------------------------PKS 2155-304------------------------------------------------------
double enef[] = {0.171, 0.4955, 1.4353, 4.1575, 12.0425, 34.882, 84.7717, 101.0379, 115.5527, 137.7252, 164.1522, 197.6806, 240.5282, 255.9009, 289.6566, 348.8196, 420.0668, 465.7592, 511.1168, 615.5136, 741.2337, 839.0094, 1120.282, 1693.1660, 3211.637};
//double obsf1[] = {1.7393, 0.3692, 0.4125,  0.7499, 1.664 , 1.821, 1.577, 1.424, 1.212, 0.975, 0.478, 0.629, 0.359, 0.214, 0.391, 0.396, 0.455, 0.402 };
//vector<double> obsf(obsf1, obsf1+18);
//double sigmaf[] = {0.9534, 0.5693, 0.726772, 1.3822599, 0.153, 0.218, 0.19, 0.208, 0.258};
//ciber
double depth[] = {0, 0, 0, 0, 3.80244e-05, 0.00222202, 0.0155311, 0.0225639, 0.0303441, 0.0454049, 0.0706268, 0.127019, 0.247873, 0.303405, 0.445013, 0.741779, 1.11812, 1.3329, 1.51827, 1.85186, 2.12762, 2.27873, 2.53497, 2.70882, 2.77346};     
//direction
double l =  17.74 * pi / 180, b =  -52.25 * pi / 180;
//datapoint numbers
double na = 25;//data points

//------------------------ICMF cluster parameter--------------------------------------------
double b2 = 2 * 2; // \mu G, in yuanq = 3, in 1406.5972 =1.  I use 2 
double kl = 0.2; //kpc^{-1}//1406 is 0.18, yuanqiang is 0.1 pi. I use 0.2 
double kh = 2 * pi ;//                  3.14,      4   pi.       2 * pi
double q = -11./3;
double eta = 1;
double rmax = 300;//kpc 
int nk = ceil(20 * (3 + log10(kh / kl)));// =log10(kh / kmin);
	  
//----------------------------MW------------------------------------------------------------------------------ 
double bx = 4.6, seita0 = 49 * pi /180, rxc = 4.8, rx = 2.9;
double h_disk = 0.4, w_disk = 0.27;
double bn = 1.4,  bs = -1.1, z0 =5.3, wh = 0.2, rn = 9.22, rs =17 ;//rs>16.7
double r0 = 8.5;
//nw: MW domains
int nw = 400;

//precedure
int n2 = 5;//5*20

double fq(double q1, double kn){
	if (kn < kl) kn = kl;
	double fq1 = (pow(kh, q1) - pow(kn, q1)) / q1;
	return fq1;
}
double Fq(double kn){
	double a1 = q + 2; double b1 = q + 3; 
	double Fq1 = (fq( a1, kn) + pow(kn, 2) * fq(q, kn)) / fq( b1, kl); 	
	return Fq1;
}
double epsilon(double kn){
	double eps1 = pi * b2 * Fq(kn) / 4;
  //  cout<< eps1 <<endl;
    return eps1;
} 

double  bx3(double x3, vector<double> bran1, vector<double> bran2){
double bx31 = 0;     
for (int n = 1; n < nk + 0.5; n++ ){
double n_1 = n;
double kn = kl * 1e-3 * pow( kh / (1e-3 * kl), n_1 / nk);
double dkn = kn -  kl * 1e-3 * pow( kh / (1e-3 * kl), (n_1 - 1) / nk);
    bx31 += pow(2 * epsilon(kn) * dkn * log(1. / bran1[n] ) / pi, 1./2 ) * cos(kn * x3 + 2 * pi * bran2[n]);
//cout<<"bx31= "<< bx31 <<endl;
}
return bx31;
}
double nr (double s){
    double n0 = 0.01;
    //double rc = 100;//kpc
    double rc = 100;//kpc// beta = 2./3;
    return n0 * pow(1 + s / rc, -1);
}
double br(const double x3, vector<double> bran1, vector<double> bran2, vector<double> bran3, vector<double> bran4){
    double bx3_1 = bx3(x3, bran1, bran2);
    double bx3_2 = bx3(x3, bran3, bran4);
    return pow((pow(bx3_1, 2) + pow(bx3_2, 2)), 1./2);
}

double seita ( double s2) {//seita(s)
    double r2 = sqrt(s2 * s2 * cos(b) * cos(b) + r0 * r0 - 2 * r0 * s2 * cos(b) * cos(l));
    if (0 < l < pi / 2 || - pi / 2 < l < 0 ) {
    if (s2 * cos(b) * cos(l) > r0) {
        return pi - asin(s2 * cos(b) * sin(l) / r2);
    }
    else {
        return asin(s2 * cos(b) *  sin(l) / r2);
        }
    }
    else{
        return asin(s2 * cos(b) *  sin(l) / r2); 
        }
    }
double b_disk(double r2, double z2){//MW_disk
    double bdisk1;
    if (r2 < 3) 
        bdisk1 = 0 ;
    else if (r2 >= 3 && r2 <5) 
        bdisk1 = 0.1 ;
    else  
        bdisk1 = 5 * 0.1 / r2;
    double bdisk = bdisk1 * (1 - pow(1 + exp(-2 * (abs(z2) - h_disk) /w_disk ), -1));
        return bdisk;
                                   }
double b_diskr (double r2, double z2){
    if (r2 >= 3 && r2 <5 ) return 0;
    else	return b_disk(r2, z2) * sin(11.5 * pi /180);//maybe there is -11.5, but this is no important, so we do not need to consider it.
}
double b_diska (double r2, double z2){
    if (r2 >=3 && r2 <5)  return  b_disk(r2, z2) * cos(0);
    else	return b_disk(r2, z2) * cos(11.5 * pi /180);
}

double b_halo(double r2, double z2){	
    if (z2 > 0) 
    return exp(-abs(z2) / z0 ) * pow(1 +exp(-2 * (abs(z2) - h_disk) /w_disk ), -1)  * bn * (1 - pow(1 + exp(-2 * (abs(r2) - rn) /wh ), -1));
    else
    return exp(-abs(z2) / z0 ) * pow(1 +exp(-2 * (abs(z2) - h_disk) /w_disk ), -1) * bs * (1 - pow(1 + exp(-2 * (abs(r2) - rs) /wh ), -1));//azimuthal = b_haloa
}

double oop(double r2, double z2){
    if(r2 > 4.8)	{
        double rp = r2 - abs(z2) / tan(seita0);
        if(rp > rxc){
            double bxrp = bx * exp(- rp/rx);  
            double oop1 = bxrp * rp / r2; 
            //cout<<"rp= "<<rp<<" oop1= "<<oop1<<endl;
            return oop1;
        }
        else return 0;
    }
    else {
        double rp = r2 * rxc / (rxc + abs(z2)/ tan(seita0));
        double bxrp = bx * exp(- rp/rx);  
        double oop1 = bxrp * pow(rp /r2, 2 );
        return oop1;
    }
}
double oopr(double r2, double z2){
    if(z2>0){   
        if (r2 > 4.8){
            return oop(r2, z2) * cos(seita0);
        }
        else {
            double rp = r2 * rxc / (rxc + abs(z2)/ tan(seita0));
            double seitax = atan(abs(z2) / (r2 - rp));
            return oop(r2 ,z2) * cos(seitax);
        }
    }
    else{
        if (r2 > 4.8){
            return -oop(r2, z2) * cos(seita0);
        }
        else {
            double rp = r2 * rxc / (rxc + abs(z2)/ tan(seita0));
            double seitax = atan(abs(z2) / (r2 - rp));
            return -oop(r2 ,z2) * cos(seitax);
        }
    }
}
double b_z(double r2, double z2){
    if (r2 > 4.8){
        return oop(r2, z2) * sin(seita0);
    }
    else {
        double rp = r2 * rxc / (rxc + abs(z2)/ tan(seita0));
        double seitax = atan(abs(z2) / (r2 - rp));
        return oop(r2 ,z2) * sin(seitax);
    }
}// only z!!!
double b_r(double r2, double z2) {
    return b_diskr(r2 ,z2) + oopr(r2 , z2);
}// disk + oopr
double b_a(double r2, double z2) {
    return b_diska(r2 ,z2) + b_halo(r2 , z2);
}
double b_x (double r2, double z2) {// bx(s) 
                                   double s2 = z2 / sin(b);
                                   return b_r(r2, z2) * sin(seita(s2)) + b_a(r2, z2) * cos (seita(s2));
                                  }

double b_y (double r2, double z2) {// by(s)
    double s2 = z2 / sin(b);
    return b_r(r2, z2) * cos(seita(s2)) - b_a(r2, z2) * sin (seita(s2));
    }
double bpre (double r2, double z2) {
    return sqrt(pow(b_x(r2, z2), 2) + pow(b_y(r2, z2), 2) + pow(b_z(r2, z2), 2) - pow((b_x(r2, z2) * sin(l) * cos(b) - b_y(r2, z2) * cos(l) * cos(b) + b_z(r2, z2) * sin(b)) ,2));
}//r_core = 10 kpc , r_max = 300 kpc;
double psi(double r2, double z2, double smax){
    double b31 = bpre(r2, z2);
    double s2_pri = 1. / 2 * smax / nw; //domain
    double z2_pri = s2_pri * sin(b);
    double r2_pri = sqrt(s2_pri * s2_pri * cos(b) * cos(b) + r0 * r0 - 2 * r0 * s2_pri * cos(b) * cos(l));
    double b32 = bpre(r2_pri, z2_pri);
    double b33 = sqrt(pow(b_x(r2, z2) - b_x(r2_pri, z2_pri), 2) + pow(b_y(r2, z2) - b_y(r2_pri, z2_pri), 2) + pow(b_z(r2, z2) - b_z(r2_pri, z2_pri), 2) - pow((b_x(r2, z2) -b_x(r2_pri, z2_pri)) * sin(l) * cos(b) - (b_y(r2, z2) - b_y(r2_pri, z2_pri)) * cos (l) * cos(b) + (b_z(r2, z2) - b_z(r2_pri, z2_pri)) * sin(b), 2) );
    return acos((pow(b31, 2) + pow(b32, 2) - pow(b33, 2)) / (2 * b31 *b32));      
}
//----------------------------------------pa-----------------------------------
int ns = kh * rmax;
double pa(double g_agamma, double ma, int i, vector<double> fai, vector<vector<double>> bran1, vector<vector<double>> bran2, vector<vector<double>> bran3, vector<vector<double>> bran4){
    double fai0 = 0;
    double s3 = 1./2 * rmax / ns; 
    double x3 = rmax / ns;
    double n1 = nr(s3) * 1e6 * 1. / (5.05e12) * 1. / (5.05e15) * 3e19;// 1/ Mev * 1 / Gev * kpc  units: kpc^{-1}
    double t1 = - 4 * pi * 1./ 137 * n1 / enef[i];//u1 = delta_pl
    double u1 = t1;  
	vector<double> bran11 = bran1[0];
	vector<double> bran22 = bran2[0];
	vector<double> bran33 = bran3[0];
	vector<double> bran44 = bran4[0];
    double rc = 100;
    double v1 = g_agamma / 1e-11 * br(s3, bran11, bran22, bran33, bran44) * 1.52e-2 * pow(1 + s3 / rc, -eta );
    double s1 = - pow(ma / 1e-9, 2) /  enef[i] * 7.8e-2  ;//kpc^{-1}
    double d1 = sqrt(pow(s1 - u1, 2) + 4 * pow(v1, 2));
    double lambda4 = t1;
    double lambda5 = 1. / 2 * (s1 + u1 - d1);
    double lambda6 = 1. / 2 * (s1 + u1 + d1); 
    double d2 = (s1-u1 +d1)/(2*d1), d3 = v1/d1, d4 = (-s1+u1 +d1)/(2*d1), d5 = lambda4 * x3, d6 = lambda5 * x3, d7 = lambda6 * x3;
    vector <complex<double>> t111 ={
        polar(1., d5), 0, 0, 0, 0, 0, 0, 0, 0};
    matrix T111(3, 3, t111);
    vector <complex<double>> t112 ={
        0, 0, 0, 0, polar(d2, d6), polar(-d3, d6), 0, polar(-d3, d6), polar(d4, d6)};
    matrix T112(3, 3, t112);
    vector <complex<double>> t113 ={
        0, 0, 0, 0, polar(d4, d7), polar(d3, d7), 0, polar(d3, d7), polar(d2, d7)};
    matrix T113(3, 3, t113);
    matrix T11 = T111 + T112 + T113; 
    for (int j = 1; j < ns ; j++){//we adopt j = ns
    s3= (j + 1. /2) * 300 / ns; //here use s3, function use s.
    n1 = nr(s3) *1e6 * 1 / (5.05e12) * 1 / (5.05e15) * 3e19;// 1/ Mev * 1 / Gev * kpc  units: kpc^{-1}
    t1 = - 4 * pi * 1./ 137 * n1 / enef[i];
    u1 = t1;
	 vector<double> bran11 = bran1[j];
	 vector<double> bran22 = bran2[j];
	 vector<double> bran33 = bran3[j];
	 vector<double> bran44 = bran4[j];
    v1 = g_agamma / 1e-11 * br(s3, bran11, bran22, bran33, bran44) * 1.52e-2 * pow(1 + s3 /rc, -1);
//    cout<<"j= "<<j<<" "<<br(s3, bran11, bran22, bran33, bran44)<<endl;
    d1 = sqrt(pow(s1 - u1, 2) + 4 * pow(v1, 2));
    lambda4 = t1;
    lambda5 = 1. / 2 * (s1 + u1 - d1);
    lambda6 = 1. / 2 * (s1 + u1 + d1);
    d2 = (s1-u1 +d1)/(2*d1), d3 = v1/d1, d4 = (-s1+u1 +d1)/(2*d1), d5 = lambda4 * x3, d6 = lambda5 * x3, d7 = lambda6 * x3;
    vector <complex<double>> v ={
        cos(fai[j]), sin(fai[j]), 0, -sin(fai[j]), cos(fai[j]), 0, 0, 0, 1};
    matrix V(3, 3, v);                        
    vector <complex<double>> t121 ={
        polar(1., d5), 0, 0, 0, 0, 0, 0, 0, 0};
    matrix T121(3, 3, t121);
    T121 = V * T121 *V.conjugate();
    vector <complex<double>> t122 ={
        0, 0, 0, 0, polar(d2, d6), polar(-d3, d6), 0, polar(-d3, d6), polar(d4, d6)};
    matrix T122(3, 3, t122);
    T122 = V * T122 *V.conjugate();
    vector <complex<double>> t123 ={
        0, 0, 0, 0, polar(d4, d7), polar(d3, d7), 0, polar(d3, d7), polar(d2, d7)};
    matrix T123(3, 3, t123);
    T123 = V * T123 *V.conjugate();
    matrix T12 = T121 + T122 + T123; 
    T11 =T12 * T11;
    }
//------------------------------------------------EBL-----------------------------------------------------
double tau1 = exp(- depth[i]/2);//new change!!!
vector <complex<double>> t31 ={
    tau1, 0, 0, 0, tau1, 0, 0, 0, 1
};
matrix T31(3, 3, t31);
//------------------------------------------------MW------------------------------------------------------
double rmax = 20;//kpc. different with hmr cal.
double smax = (sqrt(rmax * rmax - r0 * r0 * sin(l) * sin(l)) + r0 * cos(l)) / cos (b);
double s2 = 1. /2 * smax / nw; //domain
double x4 =  smax / nw; //domain
double z2 = s2 * sin(b);
double r2 = sqrt(s2 * s2 * cos(b) * cos(b) + r0 * r0 - 2 * r0 * s2 * cos(b) * cos(l));
double n = 0.02e6 * 1. / (5.05e12) * 1. / (5.05e15) * 3e19;// 1/ Mev * 1 / Gev * kpc  units: kpc^{-1} in the MW, we consider the electron density is constant, since it does not important.
double t = - 4 * pi * 1./ 137 * n / enef[i];//u = delta_pl
double u = t;  
double v = g_agamma / 1e-11 * bpre(r2, z2) * 1.52e-2;
double s = - pow(ma / 1e-9, 2) /  enef[i] * 7.8e-2  ;//kpc^{-1}
double d = sqrt(pow(s - u, 2) + 4 * pow(v, 2));
double lambda1 = t;
double lambda2 = 1. / 2 * (s + u - d);
double lambda3 = 1. / 2 * (s + u + d); 
double p = 0;
double d10 = (s-u +d)/(2*d), d11 = v/d, d12 = (-s+u +d)/(2*d), d13= lambda1 * x4, d14 = lambda2 * x4, d15 = lambda3 * x4;
vector <complex<double>> t211 ={
    polar(1., d13), 0, 0, 0, 0, 0, 0, 0, 0};
matrix T211(3, 3, t111);
vector <complex<double>> t212 ={
    0, 0, 0, 0, polar(d10, d14), polar(-d11, d14), 0, polar(-d11, d14), polar(d12, d14)};
matrix T212(3, 3, t212);
vector <complex<double>> t213 ={
    0, 0, 0, 0, polar(d12, d15), polar(d11, d15), 0, polar(d11, d15), polar(d10, d15)};
matrix T213(3, 3, t213);
matrix T21 = T211 + T212 + T213; 
for (int ms = 1; ms < nw ; ms++){//domain
    s2 = (ms + 1. /2) * smax / nw; //domain
    z2 = s2 * sin(b);
    r2 = sqrt(s2 * s2 * cos(b) * cos(b) + r0 * r0 - 2 * r0 * s2 * cos(b) * cos(l));
    v = g_agamma / 1e-11 * bpre(r2, z2) * 1.52e-2;
    d = sqrt(pow(s - u, 2) + 4 * pow(v, 2));
    lambda1 = t;
    lambda2 = 1. / 2 * (s + u - d);
    lambda3 = 1. / 2 * (s + u + d); 
    p = psi(r2, z2, smax);
    d10 = (s-u +d)/(2*d), d11 = v/d, d12 = (-s+u +d)/(2*d), d13= lambda1 * x4, d14 = lambda2 * x4, d15 = lambda3 * x4;
    vector <complex<double>> v2 ={
        cos(p), sin(p), 0, -sin(p), cos(p), 0, 0, 0, 1};
    matrix V2(3, 3, v2);                        
    vector <complex<double>> t221 ={
        polar(1., d13), 0, 0, 0, 0, 0, 0, 0, 0};
    matrix T221(3, 3, t221);
    T221 = V2 * T221 *V2.conjugate();
    vector <complex<double>> t222 ={
        0, 0, 0, 0, polar(d10, d14), polar(-d11, d14), 0, polar(-d11, d14), polar(d12, d14)};
    matrix T222(3, 3, t222);
    T222 = V2 * T222 *V2.conjugate();
    vector <complex<double>> t223 ={
        0, 0, 0, 0, polar(d12, d15), polar(d11, d15), 0, polar(d11, d15), polar(d10, d15)};
    matrix T223(3, 3, t223);
    T223 = V2 * T223 * V2.conjugate();
    matrix T22 = T221 + T222 + T223; 
    T21 =T22 * T21;
} 
T11 =T21 * T31 * T11;//whole transfer matrix
vector <complex<double>> rho={
    1./2, 0, 0, 0, 1./2, 0, 0, 0, 0
};
matrix Rho(3, 3, rho);
T11 = T11 * Rho * T11.conjugate();
//cout<<"p_sur= "<<  real(T11[0][0] + T11[1][1])<<endl;
return real(T11[0][0] + T11[1][1]);
} 
int main()
{
    random_device rd;  //Will be used to obtain a seed for the random number engine
    //n2 MF scenarios
    vector<vector<double>> fai;
    fai.resize(n2);
    for (int k1 = 0; k1 < n2; k1++){
        fai[k1].resize(ns);//每行为ns列 domain of cluster MF.
    }
    for(int j3 = 0; j3 < n2; j3++){
        for(int j=0; j < ns; j++){
            srand(time(NULL));
            fai[j3][j] = ((double)rand() / RAND_MAX) * pi;
//            cout<<"fai"<<j3<<" "<<j<< fai[j3][j]<<endl;
        }
    }
//    cout<<"fai= "<<fai[0][0]<<endl;
vector<vector<vector<double>>> bran11, bran22, bran33, bran44;
    bran11.resize(n2);
    bran22.resize(n2);
    bran33.resize(n2);
    bran44.resize(n2);
    for (int k1 = 0; k1 < n2; k1++){
        bran11[k1].resize(ns);
        bran22[k1].resize(ns);
        bran33[k1].resize(ns);
        bran44[k1].resize(ns);
    }
	 for (int k1 = 0; k1 < n2; k1++){
    for (int k2 = 0; k2 < ns; k2++){
        bran11[k1][k2].resize(nk+1);
        bran22[k1][k2].resize(nk+1);
        bran33[k1][k2].resize(nk+1);
        bran44[k1][k2].resize(nk+1);
    }
	 }
    for(int j3 = 0; j3 < n2; j3++){
        for(int j = 0; j < ns; j++){
		  for(int k = 0; k < nk +1; k++){
            srand(time(NULL));
            bran11[j3][j][k] = (double) rand() / RAND_MAX;  
            bran22[j3][j][k] = (double) rand() / RAND_MAX;  
            bran33[j3][j][k] = (double) rand() / RAND_MAX;  
            bran44[j3][j][k] = (double) rand() / RAND_MAX;  
        }
     }
	 }
//    cout<<" bran11= "<<bran11[0][0][0]<<endl;
//   cout<<" bran22= "<<bran22[0][0][0]<<endl;
double p1; 
    for(double log_ma = - 10; log_ma < -5.91; log_ma +=0.1){                       
        for(double log_gag = -12; log_gag < -8.91; log_gag +=0.1){
        double g_agamma = pow(10, log_gag), ma = pow(10, log_ma);
            for(int j4 = 0; j4 < n2; j4++){//pp scenario	
            vector<double> fai1 =  fai[j4];
            vector<vector<double>> bran1 = bran11[j4]; 
            vector<vector<double>> bran2 = bran22[j4]; 
            vector<vector<double>> bran3 = bran33[j4]; 
            vector<vector<double>> bran4 = bran44[j4]; 
            for(int i1= 0; i1 < na; i1++){
                 p1 = pa(g_agamma, ma, i1, fai1, bran1, bran2, bran3, bran4);
                 cout <<p1<<" ";
//                    pp[j7][j4][i1] = 1;
//                    cout<<j7<<" "<<j4<<" "<<i1<<endl;
                }
                    cout<<endl;
            }
        }
    }//devide job into several indpendent work, 
    }
