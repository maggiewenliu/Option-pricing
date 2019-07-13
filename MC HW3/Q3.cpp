//
//  main.cpp
//  monte carlo 3(2)
//
//  Created by Maggie on 16/2/27.
//  Copyright © 2016年 Maggie. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;

#define a0 2.50662823884
#define a1 -18.61500062529
#define a2 41.39119773534
#define a3 -25.44106049637
#define b0 -8.47351093090 
#define b1 23.08336743743
#define b2 -21.06224101826
#define b3 3.13082909833
#define c0 0.3374754822726147 
#define c1 0.9761690190917186
#define c2 0.1607979714918209
#define c3 0.0276438810333863
#define c4 0.0038405729373609
#define c5 0.0003951896511919 
#define c6 0.0000321767881768
#define c7 0.0000002888167364
#define c8 0.0000003960315187

#define min(a, b)  (((a) < (b)) ? (a) : (b))
#define max(a, b)  (((a) > (b)) ? (a) : (b))

double CDF(double x);

double r = 0.005;
double K = 100;
double S = 100;
double T = 0.01923;
double Vol = 0.3;
double no_of_trials = 20000000;
double put = 1.65466;
double Sim_p, mean, y,abs_err,std_err;


double get_uniform()
{
    return (((float) random())/(pow(2.0, 31.0)-1.0));
}


double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
}

double accep_rej()
{
    double U1 = get_uniform();
    double U2 = get_uniform();
    double x = -log(U1);
    while (U2>exp(-0.5*(x-1)*(x-1)))
    {
        U1 = get_uniform();
        U2 = get_uniform();
        x = -log(U1);
    }
    double U3 = get_uniform();
    if (U3 <= 0.5)
    {
        return -x;
    }
    else
    {
        return x;
    }
        
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

double inverse_N()
{
    double U=get_uniform();
    double r,x;
    double y=U-0.5;
    if (abs(y)<0.42)
    {
        r=y*y;
        x=y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
    }
    else
    {   r=U;
        if(y>0) r=1-U;
        r=log(-log(r));
        x=c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
        if(y<0) x=-x;
    }
    x=x-(CDF(x)-U)*exp(-x*x/2+0.5*log(2*3.14159265357));
    return x;
    
}



void box_muller()
{
    clock_t start_time=clock();
    mean = 0.0;
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = get_gaussian();
        Sim_p =S*exp((r-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T));
        double Sim_Put= exp(-r*T)*(max(Sim_p-K,0));
        mean += Sim_Put;
        y += Sim_Put*Sim_Put;
    }
    mean = mean/no_of_trials;
    y = y/no_of_trials;
    std_err=sqrt((1/(no_of_trials-1))*(y-mean*mean));
    abs_err = abs(mean-put);
    clock_t end_time=clock();
    cout << "box muller method"<<endl;
    cout << "estimate = "<< mean << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
    
}


void acceptance_rejection()
{
    clock_t start_time=clock();
    mean = 0.0;
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = accep_rej();
        Sim_p =S*exp((r-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T));
        double Sim_Put= exp(-r*T)*(max(Sim_p-K,0));
        mean += Sim_Put;
        y += Sim_Put*Sim_Put;
    }
    mean = mean/no_of_trials;
    y = y/no_of_trials;
    std_err=sqrt((1/(no_of_trials-1))*(y-mean*mean));
    abs_err = abs(mean-put);
    clock_t end_time=clock();
    cout << "acceptance_rejection method"<<endl;
    cout << "estimate = "<< mean << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
}

void inverse_transform()
{
    clock_t start_time=clock();
    mean = 0.0;
    y = 0.0;
    for(int i =0;i<no_of_trials;i++)
    {
        double ran = inverse_N();
        Sim_p =S*exp((r-0.5*pow(Vol,2))*T+Vol*ran*sqrt(T));
        double Sim_Put= exp(-r*T)*(max(Sim_p-K,0));
        mean += Sim_Put;
        y += Sim_Put*Sim_Put;
    }
    mean = mean/no_of_trials;
    y = y/no_of_trials;
    std_err=sqrt((1/(no_of_trials-1))*(y-mean*mean));
    abs_err = abs(mean-put);
    clock_t end_time=clock();
    cout << "inverse_transform method"<<endl;
    cout << "estimate = "<< mean << endl << "absolute error = "<<abs_err<<endl<<"standard error = "<<std_err<<endl;
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    cout<< "total computation time is: "<<time<<" second"<<endl;
}

int main(int argc, const char * argv[])
{
    cout << "Expiration Time (Years) = " << T << endl;
    cout << "Risk Free Interest Rate = " << r << endl;
    cout << "Volatility (%age of stock value) = " << Vol*100 << endl;
    cout << "Initial Stock Price = " << S << endl;
    cout << "Strike Price = " << K << endl;
    cout << "Price according Black Scholes  = " << put << endl;
    cout << "Number of trials : "<<no_of_trials<<endl;
    cout << "--------------------------------------" << endl;
    
    acceptance_rejection();
    cout << "--------------------------------------" << endl;
    inverse_transform();
    cout << "--------------------------------------" << endl;
    box_muller();
}
