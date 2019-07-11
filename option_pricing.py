
import math

def american_call_option(k, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob):
    if k == no_division: 
        return max(0.0, (current_stock_price - strike_price))

    up_american_call = american_call_option(k+1, i+1 , current_stock_price*up_factor, strike_price, up_prob, down_prob, notick_prob)
    down_american_call = american_call_option(k+1, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob)
    notick_american_call = american_call_option(k+1, i-1, current_stock_price/up_factor, strike_price, up_prob, down_prob, notick_prob)

    return max( (current_stock_price - strike_price), (uptick_prob * up_american_call + notick_prob*notick_american_call + down_prob * down_american_call ) / R )


def american_put_option(k, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob):
    if k == no_division: 
        return max(0.0, (strike_price - current_stock_price))

    up_american_call = american_put_option(k+1, i+1 , current_stock_price*up_factor, strike_price, up_prob, down_prob, notick_prob)
    down_american_call = american_put_option(k+1, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob)
    notick_american_call = american_put_option(k+1, i-1, current_stock_price/up_factor, strike_price, up_prob, down_prob, notick_prob)

    return max( (strike_price - current_stock_price), (uptick_prob * up_american_call + notick_prob*notick_american_call + down_prob * down_american_call ) / R )



def american_call_option_faster(k, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob, data):
    if (k, i) in data:
        return data[k,i]
    if k == no_division: 
        return max(0.0, (current_stock_price - strike_price))
    
    data[k+1, i+1] = american_call_option_faster(k+1, i+1 , current_stock_price*up_factor, strike_price, up_prob, down_prob, notick_prob, data)
    data[k+1, i] = american_call_option_faster(k+1, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob, data)
    data[k+1, i-1] = american_call_option_faster(k+1, i-1, current_stock_price/up_factor, strike_price, up_prob, down_prob, notick_prob, data)

    data[k,i] =  max( (current_stock_price - strike_price), (uptick_prob * data[k+1, i+1] + notick_prob*data[k+1, i]  + down_prob * data[k+1, i-1] ) / R )
    return data[k,i]



def american_put_option_faster(k, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob, data):
    if (k, i) in data:
        return data[k,i]

    if k == no_division: 
        return max(0.0, (strike_price - current_stock_price))

    data[k+1, i+1] = american_put_option_faster(k+1, i+1 , current_stock_price*up_factor, strike_price, up_prob, down_prob, notick_prob, data)
    data[k+1, i] = american_put_option_faster(k+1, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob, data)
    data[k+1, i-1] = american_put_option_faster(k+1, i-1, current_stock_price/up_factor, strike_price, up_prob, down_prob, notick_prob, data)

    data[k,i] =  max( (strike_price - current_stock_price ), (uptick_prob * data[k+1, i+1] + notick_prob*data[k+1, i]  + down_prob * data[k+1, i-1] ) / R )
    return data[k,i]



if __name__ == '__main__':

    expiration = 0.5
    no_division = 20
    risk_free_rate = 0.05
    vol = 0.3
    stock_price = 50
    strike_price = 40

    R = math.exp(risk_free_rate*expiration / no_division)
    up_factor = math.exp( vol * math.sqrt(2*expiration / no_division ))
    uptick_prob = ( ( R**0.5 - 1/ up_factor**0.5 ) /  (up_factor**0.5 - 1 /  up_factor**0.5 ) )**2
    downtick_prob = ( ( up_factor**0.5 -  R**0.5 ) /  (up_factor**0.5 - 1 /  up_factor**0.5 ) )**2
    notick_prob = 1 - uptick_prob - downtick_prob

    print("R = ", R)
    print("Up Factor = ", up_factor)
    print("Uptick Probability = ", uptick_prob)
    print("Downtick Probability = ", downtick_prob)
    print("Notick Probability = ", notick_prob)

    call_option_container = dict()
    call_option = american_call_option_faster(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob,call_option_container);
    print("call option ",call_option)

    put_option_container = dict()
    put_option = american_put_option_faster(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob,put_option_container);
    print("put option ",put_option)

    #call_option = american_call_option(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob);
    #print("call option ",call_option)

    #put_option = american_put_option(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob);
    #print("call option ",call_option)
