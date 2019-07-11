
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

    up_american_call = american_call_option(k+1, i+1 , current_stock_price*up_factor, strike_price, up_prob, down_prob, notick_prob)
    down_american_call = american_call_option(k+1, i, current_stock_price, strike_price, up_prob, down_prob, notick_prob)
    notick_american_call = american_call_option(k+1, i-1, current_stock_price/up_factor, strike_price, up_prob, down_prob, notick_prob)

    return max( (current_stock_price - strike_price), (uptick_prob * up_american_call + notick_prob*notick_american_call + down_prob * down_american_call ) / R )




if __name__ == '__main__':

    expiration = 0.5
    no_division = 20
    risk_free_rate = 0.05
    vol = 0.3
    stock_price = 50
    strike_price = 40

    R = math.exp(risk_free_rate*expiration / no_division)
    up_factor = math.exp( vol * math.sqrt(2*expiration / no_division ))
    uptick_prob = ( ( R**0.5 - 1/ up_factor**0.5 ) /  (up_factor**0.5 - 1 /  up_factor**0.5 ) )**0.5
    downtick_prob = ( ( up_factor**0.5 - 1/ R**0.5 ) /  (up_factor**0.5 - 1 /  up_factor**0.5 ) )**0.5
    notick_prob = 1 - uptick_prob - downtick_prob

    call_option = american_call_option(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob);
    print("call option ",call_option)

    put_option = american_put_option(0,0,stock_price,  strike_price, uptick_prob, downtick_prob, notick_prob);
    print("call option ",call_option)
