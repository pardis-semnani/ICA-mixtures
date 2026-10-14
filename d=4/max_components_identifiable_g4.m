% Finds the maximum number of components m_max that keeps the Jacobian of
% Q_{n,m,4} full-column-rank. m_proof is the bound given in the proof.
function [max_m, m_proof] = max_components_identifiable_g4(n)


Gamma_4 = floor((floor(n/2)+3) * (floor(n/2)+2)/12);
m_proof = max(Gamma_4, 1);


k = 1; %number of components
J = zeros(n^4+n^3+n^2+n,1);

while true
    A = rand(n,n);
    b = rand(n,1);
    D = rand(1,n);
    F = rand(1,n);
    if k == 1
        b1 = b;
    end

    J1 = jacobian_g4(b,A,D,F,b1);
    J = [J J1];



    if rank(J) < k*(size(J1,2))-1
        max_m = k-1;
        return
    end


    % % display(k);
    k = k+1;
end
end




