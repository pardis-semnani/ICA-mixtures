
alphaA0 = (t,r,i,a,b) -> (
    if i <= r then (
        return binomial(t-i,sub(r-i,ZZ))*b^(sub(r-i,ZZ))
    )
    else (
        return 0
    )
)

alphaB0 = (t,r,i,a,b) -> (
    if i>= r then (
        return binomial(i,sub(r,ZZ))*a^(sub(i-r,ZZ))
    )
    else (
        return 0
    )
)

AderB0 = (t,i,a,b,F) -> (
    -- F_(0,1):= 1;
    -- F_(0,2):= 1;
    return sum(for r from 0 to i-1 list F_(r,1)*F_(t-r, 2)*(i-r)*binomial(i,r)*a^(sub(i-r-1,ZZ)))
)

BderB0 = (t,i,a,b,F) -> (
    -- F_(0,1) := 1;
    -- F_(0,2) := 1;
    if i==t then (return 0);
    return sum(for r from 0 to i+1 list F_(r,1)*F_(t-r,2)*binomial(i,r-1)*(t-i)*a^(sub(i-r+1, ZZ)))
)

AderA0 = (t,i,a,b,F) -> (
    -- F_(0,1) := 1;
    -- F_(0,2) := 1;
    if i==0 then (return 0);
    return sum(for r from i-1 to t-1 list F_(r,1)*F_(t-r,2)*i*binomial(t-i,t-r-1)*b^(sub(r+1-i, ZZ)))
)

BderA0 = (t,i,a,b,F) -> (
    -- F_(0,1) := 1;
    -- F_(0,2) := 1;
    return sum(for r from i+1 to t list F_(r,1)*F_(t-r,2)*(r-i)*binomial(t-i,r-i)*b^(sub(r-i-1,ZZ)))
)



-- We denote \mu^{(1)}_{ri} with F_(i,r), \mu^{(2)}_{ri} with G_(i,r), a_1 with a, and b_2 with b.
R = QQ[a,b,F_(1,1)..F_(7,2),G_(1,1)..G_(7,2)]
F_(0,1) = 1
F_(0,2) = 1
G_(0,1) = 1
G_(0,2) = 1

-- Jacobian corresponding to \mu^{(1)}_{ri} columns for r\in [2] and i\in [7].
M1 = map(R^35, R^14, (m,n) ->
        if m<=1 and n<=1 and even n then F_(sub(1-n/2-1,ZZ),2)*alphaB0(1,n/2+1,m,a,b) 
   else if m<=1 and n<=1 and odd n then F_(sub(1-(n+1)/2,ZZ),1)*alphaB0(1,1-(n+1)/2,m,a,b) 
   else if m>=2 and m<=4 and n<=3 and even n then F_(sub(2-n/2-1,ZZ),2)*alphaB0(2,n/2+1,m-2,a,b)
   else if m>=2 and m<=4 and n<=3 and odd n then F_(sub(2-(n+1)/2,ZZ),1)*alphaB0(2,2-(n+1)/2,m-2,a,b)
   else if m>=5 and m<=8 and n<=5 and even n then F_(sub(3-n/2-1,ZZ),2)*alphaB0(3,n/2+1,m-5,a,b)
   else if m>=5 and m<=8 and n<=5 and odd n then F_(sub(3-(n+1)/2,ZZ),1)*alphaB0(3,3-(n+1)/2,m-5,a,b)
   else if m>=9 and m<=13 and n<=7 and even n then F_(sub(4-n/2-1,ZZ),2)*alphaB0(4,n/2+1,m-9,a,b)
   else if m>=9 and m<=13 and n<=7 and odd n then F_(sub(4-(n+1)/2,ZZ),1)*alphaB0(4,4-(n+1)/2,m-9,a,b)
   else if m>=14 and m<=19 and n<=9 and even n then F_(sub(5-n/2-1,ZZ),2)*alphaB0(5,n/2+1,m-14,a,b)
   else if m>=14 and m<=19 and n<=9 and odd n then F_(sub(5-(n+1)/2,ZZ),1)*alphaB0(5,5-(n+1)/2,m-14,a,b) 
   else if m>=20 and m<=26 and n<=11 and even n then F_(sub(6-n/2-1,ZZ),2)*alphaB0(6,n/2+1,m-20,a,b) 
   else if m>=20 and m<=26 and n<=11 and odd n then F_(sub(6-(n+1)/2,ZZ),1)*alphaB0(6,6-(n+1)/2,m-20,a,b) 
   else if m>=27 and m<=34 and n<=13 and even n then F_(sub(7-n/2-1,ZZ),2)*alphaB0(7,n/2+1,m-27,a,b)
   else if m>=27 and m<=34 and n<=13 and odd n then F_(sub(7-(n+1)/2,ZZ),1)*alphaB0(7,7-(n+1)/2,m-27,a,b)
   else 0);


-- Jacobian corresponding to \mu^{(2)}_{ri} columns for r\in [2] and i\in [7].
M2 = map(R^35, R^14, (m,n) ->
        if m<=1 and n<=1 and even n then G_(sub(1-n/2-1,ZZ),2)*alphaA0(1,n/2+1,m,a,b)
   else if m<=1 and n<=1 and odd n then G_(sub(1-(n+1)/2,ZZ),1)*alphaA0(1,1-(n+1)/2,m,a,b)
   else if m>=2 and m<=4 and n<=3 and even n then G_(sub(2-n/2-1,ZZ),2)*alphaA0(2,n/2+1,m-2,a,b)
   else if m>=2 and m<=4 and n<=3 and odd n then G_(sub(2-(n+1)/2,ZZ),1)*alphaA0(2,2-(n+1)/2,m-2,a,b)
   else if m>=5 and m<=8 and n<=5 and even n then G_(sub(3-n/2-1,ZZ),2)*alphaA0(3,n/2+1,m-5,a,b)
   else if m>=5 and m<=8 and n<=5 and odd n then G_(sub(3-(n+1)/2,ZZ),1)*alphaA0(3,3-(n+1)/2,m-5,a,b)
   else if m>=9 and m<=13 and n<=7 and even n then G_(sub(4-n/2-1,ZZ),2)*alphaA0(4,n/2+1,m-9,a,b)
   else if m>=9 and m<=13 and n<=7 and odd n then G_(sub(4-(n+1)/2,ZZ),1)*alphaA0(4,4-(n+1)/2,m-9,a,b)
   else if m>=14 and m<=19 and n<=9 and even n then G_(sub(5-n/2-1,ZZ),2)*alphaA0(5,n/2+1,m-14,a,b)
   else if m>=14 and m<=19 and n<=9 and odd n then G_(sub(5-(n+1)/2,ZZ),1)*alphaA0(5,5-(n+1)/2,m-14,a,b)
   else if m>=20 and m<=26 and n<=11 and even n then G_(sub(6-n/2-1,ZZ),2)*alphaA0(6,n/2+1,m-20,a,b)
   else if m>=20 and m<=26 and n<=11 and odd n then G_(sub(6-(n+1)/2,ZZ),1)*alphaA0(6,6-(n+1)/2,m-20,a,b)
   else if m>=27 and m<=34 and n<=13 and even n then G_(sub(7-n/2-1,ZZ),2)*alphaA0(7,n/2+1,m-27,a,b)
   else if m>=27 and m<=34 and n<=13 and odd n then G_(sub(7-(n+1)/2,ZZ),1)*alphaA0(7,7-(n+1)/2,m-27,a,b)
   else 0)




FF = map(R^8, R^3, (i,j) -> if j==0 then 0 else F_(i,j))
GG = map(R^8, R^3, (i,j) -> if j==0 then 0 else G_(i,j))

-- Jacobian corresponding to the column labelled by a_1.
M3 = map(R^35, R^1, (m,n) ->
        if m<=1 then AderB0(1,m,a,b,FF)
   else if m>=2 and m<=4 then AderB0(2,m-2,a,b,FF)
   else if m>=5 and m<=8 then AderB0(3,m-5,a,b,FF)
   else if m>=9 and m<=13 then AderB0(4,m-9,a,b,FF)
   else if m>=14 and m<=19 then AderB0(5,m-14,a,b,FF)
   else if m>=20 and m<=26 then AderB0(6,m-20,a,b,FF)
   else if m>=27 and m<=34 then AderB0(7,m-27,a,b,FF)
   else 0)


-- Jacobian corresponding to the column labelled by b_1.
M4 = map(R^35, R^1, (m,n) ->
        if m<=1 then BderB0(1,m,a,b,FF)
   else if m>=2 and m<=4 then BderB0(2,m-2,a,b,FF)
   else if m>=5 and m<=8 then BderB0(3,m-5,a,b,FF)
   else if m>=9 and m<=13 then BderB0(4,m-9,a,b,FF)
   else if m>=14 and m<=19 then BderB0(5,m-14,a,b,FF)
   else if m>=20 and m<=26 then BderB0(6,m-20,a,b,FF)
   else if m>=27 and m<=34 then BderB0(7,m-27,a,b,FF)
   else 0)


-- Jacobian corresponding to the column labelled by a_2.
M5  = map(R^35, R^1, (m,n) ->
        if m<=1 then AderA0(1,m,a,b,GG)
   else if m>=2 and m<=4 then AderA0(2,m-2,a,b,GG)
   else if m>=5 and m<=8 then AderA0(3,m-5,a,b,GG)
   else if m>=9 and m<=13 then AderA0(4,m-9,a,b,GG)
   else if m>=14 and m<=19 then AderA0(5,m-14,a,b,GG)
   else if m>=20 and m<=26 then AderA0(6,m-20,a,b,GG)
   else if m>=27 and m<=34 then AderA0(7,m-27,a,b,GG)
   else 0)


-- Jacobian corresponding to the column labelled by b_2.
M6 = map(R^35, R^1, (m,n) ->
        if m<=1 then BderA0(1,m,a,b,GG)
   else if m>=2 and m<=4 then BderA0(2,m-2,a,b,GG)
   else if m>=5 and m<=8 then BderA0(3,m-5,a,b,GG)
   else if m>=9 and m<=13 then BderA0(4,m-9,a,b,GG)
   else if m>=14 and m<=19 then BderA0(5,m-14,a,b,GG)
   else if m>=20 and m<=26 then BderA0(6,m-20,a,b,GG)
   else if m>=27 and m<=34 then BderA0(7,m-27,a,b,GG)
   else 0)


testValue = (i,j,a,b) -> (
    return (sum(for l from 0 to i list F_(l,1)*F_(i-l,2)* alphaB0(1,j,l,a,b)) - sum(for l from 0 to i list G_(l,1) * G_(i-l,2) * alphaA0(i,l,j,a,b)))
    )
-- Jacobian corresponding to the column labelled by \lambda_1.
M7 = map(R^35, R^1, (m,n) ->
        if m<=1 then testValue(1,m,a,b)
   else if m>=2 and m<=4 then testValue(2,m-2,a,b)
   else if m>=5 and m<=8 then testValue(3,m-5,a,b)
   else if m>=9 and m<=13 then testValue(4,m-9,a,b)
   else if m>=14 and m<=19 then testValue(5,m-14,a,b)
   else if m>=20 and m<=26 then testValue(6,m-20,a,b)
   else if m>=27 and m<=34 then testValue(7,m-27,a,b)
   else 0)



-- NN is the block of J corresponding to rows E^1:E^7 and columns labelled by lambda_1,a_1,b_1,a_2,b_2,\mu^{(k)}_{ri} for k\in [2], i\in [7], and r\in [2].
NN = map(R^35, R^33, (i,j) -> if j<=13 then M1_(i,j) 
                            else if j<=27 then M2_(i,j-14) 
                            else if j==28 then M3_(i,0) 
                            else if j==29 then M4_(i,0)
                            else if j==30 then M5_(i,0) 
                            else if j==31 then M6_(i,0)
                            else if  j==32 then M7_(i,0) )
L = (for k in toList((1,1)..(7,2)) list F_(k)=>-1) | (for k in toList((1,1)..(7,2)) list G_(k)=>-1) | {a => -1, b=>1}
NN2 = sub(NN, L)
rank NN2



