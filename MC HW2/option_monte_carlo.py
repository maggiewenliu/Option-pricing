
import math
import numpy as np
import scipy.stats
#需要import scipy.stats, 不可以只import scipy 否则会报错

def call_black_scholes(S, K, r, sigma,time, q ):
    
    
    d1 = (math.log(S/K)+(r-q)*time)/(sigma* time**0.5)+0.5*sigma* time**0.5;
    d2 = d1-(sigma* time**0.5);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1) - K*math.exp(-r*time)*scipy.stats.norm.cdf(d2);



def Monte_Carlo_Simulation(no_of_trials, T, r, S_0, vol, K, q):
    sim_option_price_sum = sim_option_price_square_sum = 0.0

    for _ in range(no_of_trials):
        
        random_var = np.random.normal(0,1,1)

        stock_price_plus = S_0 * math.exp( (r - q - 0.5*vol**2) * T + \
            vol*T**0.5 *  random_var )

        stock_price_minus = S_0 * math.exp( (r - q - 0.5*vol**2) * T - \
            vol*T**0.5 *  random_var )    

        simulated_option_price = 0.5*math.exp(-r*T)*(max(stock_price_plus-K,0)+ (max(stock_price_minus-K,0)))

        sim_option_price_sum += simulated_option_price
        sim_option_price_square_sum += simulated_option_price*simulated_option_price

    option_price = sim_option_price_sum/no_of_trials
    option_price_square = sim_option_price_square_sum / no_of_trials
    standard_error =  ( (1 / (no_of_trials - 1)) * (option_price_square - option_price*option_price) ) ** 0.5 

    upper_bound = option_price + 1.96*standard_error
    lower_bound = option_price - 1.96*standard_error

    print("option price = {0}, standard_error = {1},  Confidence Interval = [{2}, {3}] ".format(option_price,standard_error,lower_bound, upper_bound))



def Monte_Carlo_Simulation_faster(no_of_trials, T, r, S_0, vol, K, q):
    sim_option_price_sum = sim_option_price_square_sum = 0.0

    random_var = np.random.normal(0,1,no_of_trials) 

    #stock_price 是 numpy array, size = [no_of_trials, ]
    stock_price_plus = S_0 * np.exp( (r - q - 0.5*vol**2) * T + \
            vol*T**0.5 *  random_var)
    
    stock_price_minus = S_0 * np.exp( (r - q - 0.5*vol**2) * T - \
            vol*T**0.5 *  random_var)

    #np.maximum compare element wise maximum between two arrays
    # np.maximum( np.array([1,3,10]), np.array([5,6,2])) = np.array([5,6,10])
    simulated_option_price =math.exp(-r*T)* 0.5*(np.maximum( 0, stock_price_plus - strike_price )+ np.maximum( 0, stock_price_minus - strike_price))
    #compare two array maximum number
    # Size : 
    # stock_price - strike_price  is [no_of_trials, ] - [1,] =  [no_of_trials, ] 
    #simulated_option_price is  np.maximum(0, [no_of_trials, ]) = [no_of_trials, ]

    sim_option_price_sum = np.sum(simulated_option_price)
    sim_option_price_square_sum = np.sum(simulated_option_price * simulated_option_price)

    option_price = sim_option_price_sum/no_of_trials
    option_price_square = sim_option_price_square_sum / no_of_trials
    standard_error =  ( (1 / (no_of_trials - 1)) * (option_price_square - option_price*option_price) ) ** 0.5 

    print(option_price_square,standard_error )

    upper_bound = option_price + 1.96*standard_error
    lower_bound = option_price - 1.96*standard_error

    print("option price = {0}, standard_error = {1},  Confidence Interval = [{2}, {3}] ".format(option_price,standard_error,lower_bound, upper_bound))






if __name__ == '__main__':

    expiration_time=1/52
    risk_free_rate=0.003866;
    initial_stock_price=1868.99;
    volatility=0.2979;
    strike_price=1870;
    dividend_yield=0.0232;
    no_of_trials=1000000

    #Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield);
    
    Monte_Carlo_Simulation_faster(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield);
    
    
    black_scholes_call_price = call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time, dividend_yield);
    print("Black Scholes Price = ",black_scholes_call_price)
