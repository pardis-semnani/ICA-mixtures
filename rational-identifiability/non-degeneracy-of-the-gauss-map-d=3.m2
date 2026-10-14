
loadPackage "LinearTruncations"
loadPackage "Elimination"
loadPackage "GraphicalModelsMLE";



multiIndexComp = (n, d, J, I) -> (
    M = mutableMatrix map(ZZ^(n+1), n+1, (i,j)->0);
    t = 0;
    J = toList J;
    for r from 1 to n do (
        for s from 1 to n do (
            M_(r,s) = length positions(take(J, {t, t+I_(r-1)-1}), i -> i==s) 
        );
        t = t + I_(r-1);
    );
    return M;
)



-- polynomial parameterization pi_{Sym^d(C^n)} \circ P_{n,1,d}
buildPhiMap = (n,d) -> (

    L := diagonalMultidegrees(d, n);

    SSRing = QQ[for l in L list E_l];

    LL = {};
    for i from 0 to n-1 do (for j from 0 to n-1 do (if i!=j then (LL = LL | {A_(i,j)})));
    TTRing = QQ[mu_(1,1)..mu_(n,d),LL];
    for i from 0 to n-1 do (A_(i,i) = 1);
    for i from 1 to n do (mu_(i,0) = 1);


    S = set toList(1..n);
    S = S ^** d;
    for j from 1 to d-2 do (S = S/splice);
    S = toList S;
    

    

    mapPolys = {};


    for c from 0 to length L - 1 do (
        I = toList L_c;
        g = 0;
        for JJ in S do (
            MU = multiIndexComp (n, d, JJ, I);
            f =  product(for s from 1 to n list mu_(s, sum(for i from 1 to n list MU_(i,s))));
            for i from 1 to n do (
                for j from 1 to n do (
                    f = f * A_(i-1,j-1)^(MU_(i,j));
                )
            );
            g = g + f;
        );
        mapPolys = mapPolys | {g};
    );                               
  
    return map(TTRing, SSRing, mapPolys);
)







dmax = 3;
for n from 1 to 30 do ( 

-- builds the Jacobian of P_{n,1,d} at a generic point
L = {};
R = target buildPhiMap(n,dmax);
for d from 1 to dmax do (
    phiMap = buildPhiMap(n,d);
    L1 = for i in first entries vars source phiMap list sub(phiMap(i), R);
    L = L | L1;
);
J = transpose jacobian ideal L;

-- checks if the Jacobian of P_{n,1,d} is full-column-rank at theta_0
JAtTheta0 = sub(J, (flatten apply(toList(1..n), r -> {mu_(r,1)=>1, mu_(r,2)=>0, mu_(r,3)=>1})) | apply(toList((0,0)..(n-1,n-1)) - set apply(n, i->(i,i)), i -> A_i => 0));
print n; print(numcols JAtTheta0 == rank JAtTheta0))

-- chart of Grassmannian
T = frac(R);
J = sub(J, T);
M = submatrix(J, toList(0 .. numrows J - numcols J - 1), toList(0..numcols J -1))*inverse(submatrix(J, toList(numrows J - numcols J..numrows J - 1), toList(0..numcols J -1)));

-- finds the Jacobian of chart
jacobianList = apply(flatten entries M, f-> jacobianMatrixOfRationalFunction(f));
J2 = fold((g, h) -> g | h, jacobianList);

-- checks if the Jacobian of the chart is full-column-rank at the fixed point theta_0 
rank sub(J2, (flatten apply(toList(1..n), r -> {mu_(r,1)=>1, mu_(r,2)=>0, mu_(r,3)=>1})) | apply(toList((0,0)..(n-1,n-1)) - set apply(n, i->(i,i)), i -> A_i => 0)) == numcols J2

