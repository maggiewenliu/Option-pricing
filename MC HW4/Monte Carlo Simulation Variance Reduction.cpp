//
//  main.cpp
//  MCHW4
//
//  Created by Maggie on 16/3/14.
//  Copyright © 2016年 Maggie. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <time.h>
#include <fstream>
#include <cstdlib>
using namespace std;

#define max(a, b)  (((a) > (b)) ? (a) : (b))
#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define m1 2147483647
#define m2 2145483479
#define a12 63308
#define a13 -183326
#define a21 86098
#define a23 -539608
#define q12 33921
#define q13 11714
#define q21 24919
#define q23 3976
#define r12 12979
#define r13 2883
#define r21 7417
#define r23 2071
#define Invmp1 4.656612873077393e-10
int x10, x11, x12, x20, x21, x22;
double r = 0.005;
double K = 1900;
double S = 2000;
double T = 0.08333333;
double Vol = 0.3;
double q = 0.02;
double no_of_trials = 100000;
double Sim_p, mean, y,abs_err,std_err;
vector<double>yi;
vector<double>xi;
double bi;
double CDF(double x);

int Random()
{
    int h,p12,p13,p21,p23;
    x20 =x21; x21 =x22; x22 =p21 - p23; if(x22 <0) x22 =x22 +m2;
    /* Component 1 */
    h = x10/q13; p13 = -a13*(x10-h*q13)-h*r13;
    h = x11/q12; p12 = a12*(x11-h*q12)-h*r12;
    if(p13<0) p13 = p13+m1; if(p12<0) p12 = p12+m1;
    x10 = x11; x11 = x12; x12 = p12-p13; if(x12<0) x12 = x12+m1;
    /* Component 2 */
    h = x20/q23; p23 = -a23*(x20-h*q23)-h*r23;
    h = x22/q21; p21 = a21*(x22-h*q21)-h*r21;
    if(p23<0) p23 = p23+m2; if(p21<0) p21 = p21+m2;
    /*  Combination */
    if (x12<x22) return (x12-x22+m1); else return (x12-x22);}

double get_uniform1()
{
    return (((float) random())/(pow(2.0, 31.0)-1.0));
    
}

double get_uniform()
{
    int Z;
    Z=Random();
    if(Z==0)
        Z=m1;
    return(Z*Invmp1);
}

float get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}



double A_or_no(const double& S, const double& K, const double& r, const double& sigma, const double& time, const double& q)// this function is used to compute option price by Black Scholes fomula
{
    double time_sqrt = sqrt(time);
    double d1 = (log(S/K)+(r-q)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    double d2 = d1-(sigma*time_sqrt);
    return S*exp(-q*time)*CDF(d1);
}


double CDF(double x)
{
    if (x > 15.0) { return 1.0; };
    if (x < -15.0) { return 0.0; };
    double v[16]={0,1.253314137315500,0.6556795424187985,0.4213692292880545,0.3045902987101033,0.2366523829135607,0.1928081047153158,0.1623776608968675, 0.1401041834530502,0.1231319632579329, 0.1097872825783083,0.09902859647173193, 0.09017567550106468,0.08276628650136917,0.0764757610162485,0.07106958053885211};
    double c = 0.918938533204672;
    int j = floor(min(abs(x)+0.5,14));
    int z = j;
    double h = abs(x)-z;
    double a = v[j+1];
    double b = z*a-1;
    double q =1;
    double s = a+ h *b;
    for (int i = 2;i<=24-j;i = i+2)
    {
        a = (a+z*b)/i;
        b = (b+z*a)/(i+1);
        q = q*h*h;
        s = s+q*(a+h*b);
    }
    double y = s*exp(-0.5*x*x-c);
    if(x>0)
        y = 1-y;
    return y;
}

void monte_carlo_simulation(double call)// this function is built by the calculate the mean of option price and the mean of option price square. This method is to avoid calculating error by updating mean in no_of_trials times
{
    clock_t start_time=clock();
    mean = 0.0;
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = get_gaussian();
        Sim_p =S*exp((r-q-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T));
        if (Sim_p>=K)
        {
            Sim_p = exp(-r*T)*Sim_p;
            mean += Sim_p;
            y += Sim_p*Sim_p;
        }
    }
    mean = mean/no_of_trials;
    y = y/no_of_trials;
    std_err=sqrt((1/(no_of_trials-1))*(y-mean*mean));
    abs_err = abs(mean-call);
    clock_t end_time=clock();
    cout << "monte carlo simulation"<<endl;
    cout << "estimate = "<< mean << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
}


void control(double call)// this function is built by the calculate the mean of option price and the mean of option price square. This method is to avoid calculating error by updating mean in no_of_trials times
{
    double var_x = 0.0;
    double var_y = 0.0;
    clock_t start_time=clock();
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = get_gaussian();
        Sim_p =S*exp((r-q-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T));
        xi.push_back(Sim_p);
        if(Sim_p>=K)
        {
            Sim_p = exp(-r*T)*Sim_p;
            y += Sim_p;
            yi.push_back(Sim_p);
        }
        else{
            yi.push_back(0);
        }
    }
    mean = S*exp((r-q)*T);
    y = y/no_of_trials;
    double above = 0;
    double below = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        var_x = var_x + (xi[i]-mean)*(xi[i]-mean);
        var_y = var_y + (yi[i]-y)*(yi[i]-y);
        above= above + (xi[i]-mean)*(yi[i]-y);
        below = below + (xi[i]-mean)*(xi[i]-mean);
    }
    var_y = var_y/(no_of_trials-1);
    var_x = var_x/(no_of_trials-1);
    double covar = above/(no_of_trials-1);
    double pxy = covar/(sqrt(var_x)*sqrt(var_y));
    double b_hat = above/below;
    double yb = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
    yb += (yi[i-1]+b_hat*(mean-xi[i-1]));
    }
    yb = yb/no_of_trials;
    std_err=sqrt((1-pxy*pxy)*(var_y/no_of_trials));
    abs_err = abs(yb-call);
    
    clock_t end_time=clock();
    cout << "control variate simulation"<<endl;
    cout << "estimate = "<< yb << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
}

void sampling( double call)
{
    double mu = (r-q-0.5*Vol*Vol)*T;
    double mu_hat = 0.051;
    clock_t start_time=clock();
    mean = 0.0;
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = get_gaussian();
        double x=(r-q-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T)+mu_hat-mu;
        if (x>=log(K/S))
        {
            Sim_p = S*exp(x)*exp(-(mu*mu - mu_hat*mu_hat+2*(mu_hat-mu)*x)/(2*Vol*Vol*T));
            Sim_p = exp(-r*T)*Sim_p;
            mean += Sim_p;
            y += Sim_p*Sim_p;
        }
    }
    mean = mean/no_of_trials;
    y = y/no_of_trials;
    std_err=sqrt((1/(no_of_trials-1))*(y-mean*mean));
    abs_err = abs(mean-call);
    clock_t end_time=clock();
    cout << "important sampling simulation"<<endl;
    cout << "estimate = "<< mean << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
}



int main(int argc, const char * argv[]) {
    x11 = get_uniform1()*2147483646;
    x12 = get_uniform1()*2147483646;
    x10 = get_uniform1()*2147483646;
    x20 = get_uniform1()*2145483478;
    x21 = get_uniform1()*2145483478;
    x22 = get_uniform1()*2145483478;
    double call1=A_or_no(2000,1900,0.005,0.3,0.083333333,0.02);
    double call2=A_or_no(2000,2200,0.005,0.3,0.083333333,0.02);
    cout << "option price with strike 1900 is "<<call1<<endl;
    cout << "option price with strike 2200 is "<<call2<<endl;
    monte_carlo_simulation(call1);
    control(call2);
    sampling(call2);
}





