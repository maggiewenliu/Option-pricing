
import math
import numpy as np
import scipy.stats
#需要import scipy.stats, 不可以只import scipy 否则会报错

def call_black_scholes(S, K, r, sigma,time, q ):
    
    
    d1 = (math.log(S/K)+(r-q)*time)/(sigma* time**0.5)+0.5*sigma* time**0.5;
    d2 = d1-(sigma* time**0.5);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1) - K*math.exp(-r*time)*scipy.stats.norm.cdf(d2);


def A_or_no(S, K, r, sigma, time, q):
    time_sqrt = time**0.5
    d1 = (math.log(S/K)+(r-q)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    d2 = d1-(sigma*time_sqrt);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1);



def Monte_Carlo_Simulation(no_of_trials, T, r, S_0, vol, K, q,call_price):
    sim_option_price_sum = sim_option_price_square_sum = 0.0

    for _ in range(no_of_trials):

        stock_price = S_0 * math.exp( (r - q - 0.5*vol**2) * T + vol*T**0.5 *  np.random.normal(0,1,1)  )    

        if stock_price >= K:
            simulated_option_price = math.exp(-r*T)*stock_price

            sim_option_price_sum += simulated_option_price
            sim_option_price_square_sum += simulated_option_price*simulated_option_price

    option_price = sim_option_price_sum/no_of_trials
    option_price_square = sim_option_price_square_sum / no_of_trials
    std_err =  ( (1 / (no_of_trials - 1)) * (option_price_square - option_price*option_price) ) ** 0.5 
    abs_err = abs(option_price-call_price)
    
    print("Direct Approach option price = {0}, standard_error = {1},  abs_err = {2} ".format(option_price,std_err,abs_err))


def Control_Variate(no_of_trials, T, r, S_0, vol, K, q,call_price):
    x, y = [], [] #先用list 因为numpy.array append 速度太慢，跟list 在for loop 中速度慢10倍
    for _ in range(no_of_trials):
        sim_price =  S_0 * math.exp( (r - q - 0.5*vol**2) * T + vol*T**0.5 *  np.random.normal(0,1,1) )
        x.append(sim_price) #x = np.append(x, sim_price) 这是numpy array appen的方式

        if sim_price >= K:
            sim_price = math.exp(-r*T)*sim_price
            y.append(sim_price)
        else:
            y.append(0)
    x, y = np.array(x), np.array(y)

    mean = S_0 *math.exp((r-q)*T)
    y_mean = np.average(y)

    var_x = 1/(no_of_trials-1) * np.sum( (x - mean)**2 ) #把每个数减去平均数，然后平方
    var_y = 1/(no_of_trials-1) * np.sum( (y - y_mean)**2 )
    covar = 1/(no_of_trials-1) * np.sum( (x - mean) *(y - y_mean))

    pxy = covar / (var_x * var_y)**0.5
    b_har = covar/ var_x

    yb = y + b_har * ( mean - x) #yb 是numpy array, x也是numpy_array

    option_price = np.average(yb)
    std_err = ( (1-pxy*pxy) * (var_y/ (no_of_trials -1)) ) **0.5
    abs_err = abs(option_price - call_price)

    print("Control Variate option price = {0}, standard_error = {1},  abs_err = {2} ".format(option_price,std_err,abs_err))


def Importance_Sampling(no_of_trials, T, r, S_0, vol, K, q,call_price):
    
    mu = (r-q-0.5*vol*vol)*T;
    mu_hat = 0.051;

    option_price_sum = 0.0;
    option_price_sum_square = 0.0;
    
    for _ in range(no_of_trials):

        x=(r-q- 0.5*vol**2 )*T+ vol* T**0.5 *  np.random.normal(0,1,1)+mu_hat-mu;

        if x >= math.log(K/S_0):
            simulated_option_price = S_0*math.exp(x)*math.exp(-(mu*mu - mu_hat*mu_hat+2*(mu_hat-mu)*x)/(2*vol*vol*T));
            simulated_option_price = math.exp(-r*T)*simulated_option_price

            option_price_sum += simulated_option_price
            option_price_sum_square += simulated_option_price*simulated_option_price

    option_price = option_price_sum/no_of_trials
    option_price_square = option_price_sum_square / no_of_trials
    std_err =  ( (1 / (no_of_trials - 1)) * (option_price_square - option_price*option_price) ) ** 0.5 
    abs_err = abs(option_price-call_price)
    
    print("Importance Sampling option price = {0}, standard_error = {1},  abs_err = {2} ".format(option_price,std_err,abs_err))





if __name__ == '__main__':

    expiration_time=1/12
    risk_free_rate=0.005
    initial_stock_price=2000;
    volatility=0.3
    strike_price=1900
    strike_price2=2200
    dividend_yield=0.02
    no_of_trials=1000000

    call_price1 = A_or_no(initial_stock_price,strike_price,risk_free_rate,volatility,expiration_time,dividend_yield)
    call_price2 = A_or_no(initial_stock_price,strike_price2,risk_free_rate,volatility,expiration_time,dividend_yield)

    print("Asset or Nothing at ", strike_price  ," option = ",call_price1)
    print("Asset or Nothing at ", strike_price2  ," option = ",call_price2)


    Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,call_price1)
    Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price2, dividend_yield,call_price2)


    Control_Variate(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,call_price1)
    Control_Variate(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price2, dividend_yield,call_price2)


    Importance_Sampling(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,call_price1)
    Importance_Sampling(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price2, dividend_yield,call_price2)

    
   