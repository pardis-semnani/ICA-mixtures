loadPackage "LinearTruncations"
loadPackage "Elimination"


-- Jacobian of  P_{n,1,d} at ({a_{rs} = 0}_{r \neq s \in [n]}, mu = {mu_{ri}}_{r \in [n], i \in [d]}); only the rows corresponding to E^d
jacob = (n,d,mu) -> (
    L := diagonalMultidegrees(d, n);
    R := ring mu_(1,1);
    J := mutableMatrix map(R^(length L), n^2+(d-1)*n, (i,j) -> 0);
    for c from 0 to length L - 1 do (
        I = toList L_c;
        counter = 0;
        for r from 1 to n do (
            for s from 1 to n do (
                if r != s then (
                    
                    if I_(r-1) != 0 then (
                        J_(c, counter) = I_(r-1) * product(for j in (toList(1..n) - set{r,s}) list mu_(j, I_(j-1)));
                        J_(c, counter) = J_(c, counter) * mu_(r, I_(r-1) -1) * mu_(s, I_(s-1) + 1);
                    );
                    counter = counter + 1;
                )
            )
        );

        for r from 1 to n do (
            for i from 1 to d do (

                if I_(r-1) == i then (
                    J_(c, counter) = product(for j in (toList(1..n) - set{r}) list mu_(j, I_(j-1)));
                );

                counter = counter + 1;
            )
        )

    );
    return J
)




d = 4; 
for n from 1 to 7 do (

    R = QQ[mu1_(1,1)..mu1_(n,d),alpha_(1,1)..alpha_(n^2+(d-1)*n,n^2+(d-1)*n)];
    for i from 1 to n do (mu1_(i,0) = 1);

    -- specific point in the proof
    mu2 = mutableMatrix map(QQ^(n+1), d+1, (n,i) -> if n==0 then 0 else if i==0 or i==1  or i == 3 then 1 else 0); 

    if n < 3 then (
        -- the Jacobina at the fixed point
        J2 = (matrix jacob(n,d-3,mu2) | map (R^(binomial(n+d-4,d-3)),3*n,(i,j)->0)) || (matrix jacob(n,d-2,mu2) | map (R^(binomial(n+d-3,d-2)),2*n,(i,j)->0)) || (matrix jacob(n,d-1,mu2) | map (R^(binomial(n+d-2,d-1)),n,(i,j)->0)) || matrix jacob(n,d,mu2);
        
        -- the Jacobian at a generic point
        J1 = (matrix jacob(n,d-3,mu1) | map (R^(binomial(n+d-4,d-3)),3*n,(i,j)->0)) || (matrix jacob(n,d-2,mu1) | map (R^(binomial(n+d-3,d-2)),2*n,(i,j)->0)) || (matrix jacob(n,d-1,mu1) | map (R^(binomial(n+d-2,d-1)),n,(i,j)->0)) || matrix jacob(n,d,mu1);
    )

    else (
        -- the Jacobian at the fixed point (only the rows corresponding to E^3, E^4)
        J2 = (matrix jacob(n,d-1,mu2) | map (R^(binomial(n+d-2,d-1)),n,(i,j)->0)) || matrix jacob(n,d,mu2); 

        -- the Jacobian at a generic point (only the rows corresponding to E^3, E^4)
        J1 = (matrix jacob(n,d-1,mu1) | map (R^(binomial(n+d-2,d-1)),n,(i,j)->0)) || matrix jacob(n,d,mu1);
    );
    
    
    -- check if the Jacobian of P_{n,1,d-1} at the fixed point is full-column-rank
    J3 = (matrix jacob(n,d-3,mu2) | map (R^(binomial(n+d-4,d-3)),2*n,(i,j)->0)) || (matrix jacob(n,d-2,mu2) | map (R^(binomial(n+d-3,d-2)),n,(i,j)->0)) || (matrix jacob(n,d-1,mu2));
    print(rank(J3) == numcols(J3));



    a = map(R^(n^2+(d-1)*n), R^(n^2+(d-1)*n), (i,j) -> alpha_(i+1,j+1));
    I = ideal flatten entries (J2*a - J1);
    I2 = eliminate(toList(alpha_(1,1)..alpha_(n^2+(d-1)*n,n^2+(d-1)*n)), I);

    -- check if the fixed point is recoverable from the column space of the Jacobian, up to a finite class. 
    if n*d == codim I2 then (print (n,d)) else (print codim I2);
)
      













