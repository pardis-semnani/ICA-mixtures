% Finds the maximum number of components m_max that keeps the Jacobian of
% Q_{n,m,3} full-column-rank.
function max_m = max_components_identifiable_g3(n)

k = 1; %number of components
J = zeros(n^3+n^2+n,1);
while true
    A = rand(n,n);
    D = rand(1,n);
    b = rand(n,1);
    if k==1
        b1 = b;
    end


    J = [J jacobian_g3(b,A,D,b1)];

    if rank(J) < k*(2*n+n^2+1)-1
        max_m = k-1;
        return
    end
k = k+1;
end
end