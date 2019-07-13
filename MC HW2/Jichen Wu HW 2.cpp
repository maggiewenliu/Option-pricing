//
//  main.cpp
//  j
//
//  Created by Maggie on 16/1/24.
//  Copyright © 2016年 Maggie. All rights reserved.
//

//
//  main.cpp
//  hw1 E525
//
//  Created by Maggie on 16/1/23.
//  Copyright © 2016年 Maggie. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include<time.h>
using namespace std;

double standard_error;

#define max(a, b)  (((a) > (b)) ? (a) : (b))
class call_option
{
    double risk_free_rate,strike_price,initial_stock_price, expiration_time,volatility,dividend_yield;;
    double no_of_trials;
    double price_plus,price_minus,sim_option_price_sum, option_price_square,sim_option_price_square_sum;
    double option_price,lower_bound,upper_bound;
    double no_of_trials_to_get_two_cent_CI; // this variable is compute the number of trials to get the confidence interval that is 2 cent wide
    float get_uniform()
    {
        return (((float) random())/(pow(2.0, 31.0)-1.0));
    }
    
    float get_gaussian()
    {
        return (sqrt(-2.0*log(get_uniform()))*cos(6.283185307999998*get_uniform()));
    }
    
    void monte_carlo_simulation(double no_of_trials)// this function is built by the calculate the mean of option price and the mean of option price square. This method is to avoid calculating error by updating mean in no_of_trials times
    {
        for(int i =0;i<no_of_trials;i++)
        {
            double random = get_gaussian();
            price_plus = initial_stock_price*exp((risk_free_rate-dividend_yield-0.5*pow(volatility,2))*expiration_time+volatility*random*sqrt(expiration_time));
            price_minus=initial_stock_price*exp((risk_free_rate-dividend_yield-0.5*pow(volatility,2))*expiration_time-volatility*random*sqrt(expiration_time));
            double simulated_option_price= 0.5*(max(price_plus-strike_price,0.0)+max(price_minus-strike_price,0.0))*exp(-risk_free_rate*expiration_time);
            sim_option_price_sum += simulated_option_price;
            sim_option_price_square_sum += simulated_option_price*simulated_option_price;
        }
        option_price = sim_option_price_sum/no_of_trials;
        option_price_square = sim_option_price_square_sum/no_of_trials;
        standard_error=sqrt((1/(no_of_trials-1))*(option_price_square-option_price*option_price));
        upper_bound = option_price+1.96*standard_error;
        lower_bound = option_price-1.96*standard_error;
    }
    
    void monte_carlo_simulation_2(double no_of_trials)// this function is built by the update the mean of option square and the mean of option price square in the process of simulation
    {
        for(int i =0;i<no_of_trials;i++)
        {
            int random = get_gaussian();
            price_plus=initial_stock_price*exp((risk_free_rate-dividend_yield-0.5*pow(volatility,2))*expiration_time+volatility*random*sqrt(expiration_time));
            price_minus=initial_stock_price*exp((risk_free_rate-dividend_yield-0.5*pow(volatility,2))*expiration_time-volatility*random*sqrt(expiration_time));
            double simulated_option_price= 0.5*exp(-risk_free_rate*expiration_time)*(max(price_plus-strike_price,0.0)+max(price_minus-strike_price,0.0));
            option_price = (option_price *i + simulated_option_price)/(i+1);
            option_price_square = (option_price_square *i + simulated_option_price*simulated_option_price)/(i+1);
        }
        standard_error=sqrt((1/(no_of_trials-1))*(option_price_square-option_price*option_price));
        upper_bound = option_price+1.96*standard_error;
        lower_bound = option_price-1.96*standard_error;
    }
    

    
    void print_result()
    {
        cout << "Expiration Time (Years) = " << expiration_time << endl;
        cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
        cout << "Volatility (%age of stock value) = " << volatility*100 << endl;
        cout << "Dividend Rate = "<<dividend_yield<<endl;
        cout << "Initial Stock Price = " << initial_stock_price << endl;
        cout << "Strike Price = " << strike_price << endl;
        cout << "number of trials = " << no_of_trials << endl;
        cout << "---------------------------------"<<endl;
        cout <<"Monte Carlo European Option Pricing = " <<option_price<<endl;
        cout <<"Monte Carlo Standard Error = "<<standard_error<<endl;
        cout <<"95% confidence Interval = "<< "["<<lower_bound<<","<<upper_bound<<"]"<<endl;
    }
    
public:
    void main(int argc, char * argv[])
    {
        sscanf(argv[1], "%lf",&no_of_trials);
        expiration_time=0.01923;
        risk_free_rate=0.003866;
        initial_stock_price=1868.99;
        volatility=0.2979;
        strike_price=1870;
        dividend_yield=0.0232;
        //no_of_trials=100000;
        sim_option_price_sum = 0.0;
        sim_option_price_square_sum=0.0;
        option_price = 0.0;
        option_price_square = 0.0;
        monte_carlo_simulation(no_of_trials);
        print_result();
    }
};

int main(int argc, char* argv[])
{
    clock_t start_time=clock();
    call_option x;
    x.main(argc, argv);
    clock_t end_time=clock();
    double time = static_cast<double>(end_time-start_time)/CLOCKS_PER_SEC;
    double efficiency = time * pow(standard_error,2);
    cout<< "Running time is: "<<time<<" second"<<endl;
    cout << " efficiency is "<< efficiency << endl;
    return 0;
}
















