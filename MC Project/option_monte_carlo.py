
import math
import numpy as np
import scipy.stats
#需要import scipy.stats, 不可以只import scipy 否则会报错
import sobol_seq

def call_black_scholes(S, K, r, sigma,time, q ):
    
    
    d1 = (math.log(S/K)+(r-q)*time)/(sigma* time**0.5)+0.5*sigma* time**0.5;
    d2 = d1-(sigma* time**0.5);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1) - K*math.exp(-r*time)*scipy.stats.norm.cdf(d2);


def geometric_black_scholes_call(S,  K, r, sigma, time, q): # this function is used to compute option price by Black Scholes fomula
    time_sqrt = time ** 0.5
    d1 = (math.log(S/K)+(r-q)*time)/(sigma*time_sqrt)+0.5*sigma*time_sqrt;
    d2 = d1-(sigma*time_sqrt);
    return S*math.exp(-q*time)*scipy.stats.norm.cdf(d1) - K*math.exp(-r*time)*scipy.stats.norm.cdf(d2);


def Monte_Carlo_Simulation(no_of_trials, T, r, S_0, vol, K,  q, m):#m是分几个interval
    sim_option_price_sum = sim_option_price_square_sum = 0.0
    itval = T/m;

    for _ in range(no_of_trials):

        average = price = S_0*math.exp((r-q-0.5* vol**2 )*itval+vol*np.random.normal(0,1,1)* itval **0.5 );

        for _ in range(m-1):
            price = price*math.exp((r-q-0.5* vol**2 )*itval+vol*np.random.normal(0,1,1)* itval **0.5 );
            average += price 
        
        average = average / m;
        simulated_option_price = math.exp(-r*T)*max(average-K, 0)
        sim_option_price_sum += simulated_option_price
        sim_option_price_square_sum += simulated_option_price*simulated_option_price

        
    option_price = sim_option_price_sum/no_of_trials
    option_price_square = sim_option_price_square_sum / no_of_trials
    std_err =  ( (1 / (no_of_trials - 1)) * (option_price_square - option_price*option_price) ) ** 0.5 
    
    upper_bound = option_price + 1.96*std_err
    lower_bound = option_price - 1.96*std_err

    print("Direct Approach option price = {0}, standard_error = {1},  Confidence Interval = [{2:.3f}, {3:.3f}] ".format(option_price,std_err,upper_bound,lower_bound))


def Control_Variate(no_of_trials, T, r, S_0, vol, K, q,m):
    x, y = [], [] #先用list 因为numpy.array append 速度太慢，跟list 在for loop 中速度慢10倍
    itval = T/m;

    for _ in range(no_of_trials):
        geo_avrg = average = price = S_0*math.exp((r-q-0.5* vol**2 )*itval+vol*np.random.normal(0,1,1)* itval **0.5 );
        
        for _ in range(m-1):
            price = price*math.exp((r-q-0.5* vol**2 )*itval+vol*np.random.normal(0,1,1)* itval **0.5 );
            average += price 
            geo_avrg *= price 
        
        average = average / m;
        geo_avrg = geo_avrg**(1/m)
        arith_option = max(average-K, 0)
        geo_option = max(geo_avrg-K, 0)

        x += [geo_option] #等于append的
        y.append(arith_option)
        

    x, y = np.array(x), np.array(y)
   
    T_til =(0.5*(T+itval))
    sigma_til = (2*m+1)*vol*vol/(3*m)
    q_til = q+0.5*vol*vol-0.5*sigma_til
    x_mean =  math.exp(r*T_til)*geometric_black_scholes_call(S_0,K,r, sigma_til**0.5,T_til,q_til)
    y_mean =  np.mean(y)

    var_x = 1/(no_of_trials-1) * np.sum( (x - x_mean)**2 ) #把每个数减去平均数，然后平方
    var_y = 1/(no_of_trials-1) * np.sum( (y - y_mean)**2 )
    covar = 1/(no_of_trials-1) * np.sum( (x - x_mean) *(y - y_mean))

    pxy = covar / (var_x * var_y)**0.5
    b_hat = covar/ var_x

    yb = math.exp(-r*T)*(y + b_hat * ( x_mean - x))#yb 是numpy array, x也是numpy_array

    option_price = np.average(yb)
    std_err = ( (1-pxy*pxy) * (var_y/ no_of_trials) ) **0.5
    

    upper_bound = option_price + 1.96*std_err
    lower_bound = option_price - 1.96*std_err

    print("Control Variate option price = {0}, standard_error = {1},  Confidence Interval = [{2:.3f}, {3:.3f}] ".format(option_price,std_err,upper_bound,lower_bound))







def Quasi_Monte_Carlo_Simulation(no_of_trials, T, r, S_0, vol, K,  q, m, batch):#m是分几个interval
    
    
    def quasi_get_Gaussian(u1,u2):
        output_u1 =  math.sqrt(-2.0*math.log(u1))*math.cos(6.283185307999998*u2);
        output_u2 =  math.sqrt(-2.0*math.log(u1))*math.sin(6.283185307999998*u2);
        return [output_u1, output_u2]

    
    itval = T/m;

    quasi_mean = quasi_mean_square = 0
    option_price_list = []

    for _ in range(batch):

        sobel = sobol_seq.i4_sobol_generate(m,no_of_trials)
        #生成一个shape 为  no_of_trials * m 的np.array()

        random_var = np.zeros((no_of_trials,m))
        #生一个numpy array, shape 为no_of_trials * m

        sim_option_price_sum  = 0.0

        for i in range(no_of_trials):
            #print( random_var.shape, sobel.shape,len(random_var[i]), len(sobel[:,i]))
            random_var[i] = sobel[i] + np.random.uniform(0,1,m) #把sobel 所有第i行数据加上一个 uniform variable
            random_var[i] = random_var[i] - np.floor(random_var[i]) #让 所有第i行数据 都小于0
            #np.floor([2.8,1.2]) = np.array([2，1])
            random_var[i,0:2] = quasi_get_Gaussian(*random_var[i,0:2]) 
            # * 号表示unload list 元素到variable, 也可以写成quasi_get_Gaussian(random_var[i,0],random_var[i,1] ) 

            average = price = S_0*math.exp((r-q-0.5* vol**2 )*itval+vol*random_var[i,0]* itval **0.5 );

            for j in range(m-1):

                if j % 2 == 1:
                    random_var[i,j+1:j+3] = quasi_get_Gaussian(*random_var[i,j+1:j+3])

                price = price*math.exp((r-q-0.5* vol**2 )*itval+vol*random_var[i,j+1]* itval **0.5 );
                average += price 
            
            average = average / m;
            simulated_option_price = math.exp(-r*T)*max(average-K, 0)
            sim_option_price_sum += simulated_option_price
            option_price_list.append(average)
        
        option_price = average/no_of_trials
        quasi_mean += option_price
        quasi_mean_square += option_price *option_price
   
    quasi_mean = quasi_mean / batch
    quasi_mean_square = quasi_mean_square /  batch
    option_price_list = np.array(option_price_list)

    std_err = math.sqrt( 1/ batch * np.sum( (option_price_list - quasi_mean)*2 )/(batch-1 ))
    #把list中每个option price 减去quasi_mean 再平方


    upper_bound = quasi_mean + 1.96*std_err
    lower_bound = quasi_mean - 1.96*std_err

    print("Quasi Monte Carlo option price = {0}, standard_error = {1},  Confidence Interval = [{2:.3f}, {3:.3f}] ".format(option_price,std_err,upper_bound,lower_bound))




if __name__ == '__main__':

    expiration_time=2
    risk_free_rate=0.05
    initial_stock_price=2;
    volatility=0.5
    strike_price=2
    dividend_yield=0
    no_of_trials=10000
    no_of_intervals = 30
    batch = 10

    #Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,no_of_intervals)

    #Control_Variate(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,no_of_intervals)


    Quasi_Monte_Carlo_Simulation(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,no_of_intervals,batch)

    #Control_Variate(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,call_price1)
    #Control_Variate(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price2, dividend_yield,call_price2)


    #Importance_Sampling(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price, dividend_yield,call_price1)
    #Importance_Sampling(no_of_trials, expiration_time, risk_free_rate, initial_stock_price, volatility, strike_price2, dividend_yield,call_price2)

    
   