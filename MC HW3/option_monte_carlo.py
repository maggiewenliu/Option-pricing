
import math
import numpy as np
import scipy.stats
#需要import scipy.stats, 不可以只import scipy 否则会报错

def call_black_scholes(S, K, r, sigma,time, q ):
    
    
    d1 = (math.log(S/K)+(r-q)*time)/(sigma* time**0.5)+0.5*sigma* time**0.5;
    d2 = d1-(sigma* time**0.5);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1) - K*math.exp(-r*time)*scipy.stats.norm.cdf(d2);

def acceptance_rejection():
    U1 =np.random.uniform(0,1,1) #low is 0, high is 1, size = 1
    U2 = np.random.uniform(0,1,1)
    x = -math.log(U1);
    
    while U2>math.exp(-0.5*(x-1)*(x-1)):
        U1 = np.random.uniform(0,1,1)
        U2 = np.random.uniform(0,1,1)
        x = -math.log(U1)
    
    U3 = np.random.uniform(0,1,1)
    if (U3 <= 0.5):
        return -x;
    else:
        return x;

def box_muller():
    x = np.random.uniform(0,1,1)
    y = np.random.uniform(0,1,1)
    return math.sqrt( -2.0*math.log(x) )*math.cos(6.283185307999998*y)

def inverse_transform():

    a0 = 2.50662823884
    a1 =-18.61500062529
    a2= 41.39119773534
    a3 =-25.44106049637
    b0 =-8.47351093090 
    b1= 23.08336743743
    b2 =-21.06224101826
    b3 =3.13082909833
    c0= 0.3374754822726147 
    c1 =0.9761690190917186
    c2 =0.1607979714918209
    c3 =0.0276438810333863
    c4 =0.0038405729373609
    c5 =0.0003951896511919 
    c6 = 0.0000321767881768
    c7 =0.0000002888167364
    c8 =0.0000003960315187

    U=np.random.uniform(0,1,1) 
    y=U-0.5;
    if (abs(y)<0.42):
        r=y*y;
        x=y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
    else:
        r=U;
        if y>0: r=1-U;
        r=math.log(-math.log(r));
        x=c0+r*(c1+r*(c2+r*(c3+r*(c4+r*(c5+r*(c6+r*(c7+r*c8)))))));
        if y<0: x=-x;
    
    x=x-(scipy.stats.norm.cdf(x)-U)*math.exp(-x*x/2+0.5*math.log(2*3.14159265357));
    return x;




def Monte_Carlo_Simulation(no_of_trials, T, r, S_0, vol, K, q,method):
    print("Method ",method.__name__) #打印function 的name
    sim_option_price_sum = sim_option_price_square_sum = 0.0

    for _ in range(no_of_trials):
        
        random_var = method()
        #print(random_var)

        stock_price = S_0 * math.exp( (r - q - 0.5*vol**2) * T + vol*T**0.5 *  random_var )    

        simulated_option_price = math.exp(-r*T)*max(stock_price-K,0)

        sim_option_price_sum += simulated_option_price
        sim_option_price_square_sum += simulated_option_price*simulated_option_price

    option_price = sim_option_price_sum/no_of_trials
    option_price_square = sim_option_price_square_sum / no_of_trials
    standard_error = (1 / (no_of_trials - 1)) ** 0.5 * (option_price_square - option_price*option_price)

    upper_bound = option_price + 1.96*standard_error
    lower_bound = option_price - 1.96*standard_error

    print("option price = {0}, standard_error = {1},  Confidence Interval = [{2}, {3}] ".format(option_price,standard_error,lower_bound, upper_bound))







if __name__ == '__main__':

    expiration_time=1/52
    risk_free_rate=0.005
    initial_stock_price=100;
    volatility=0.3;
    strike_price=100
    dividend_yield=0;
    no_of_trials=1000000

    #Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield);
    
    #Monte_Carlo_Simulation_faster(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield);
    
    Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,acceptance_rejection);
    Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,inverse_transform);
    Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,box_muller);
    
    
    black_scholes_call_price = call_black_scholes(initial_stock_price, strike_price, risk_free_rate, volatility, expiration_time, dividend_yield);
    print("Black Scholes Price = ",black_scholes_call_price)
