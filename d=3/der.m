function B = der(A) %der(A(:,1),...,A(:,d))
d = size(A,2);
n = size(A,1);
P = perms(1:d+1);
B = zeros(n^(d+1),n);
for i=1:factorial(d+1)
    C = 1;
    for j=1:d+1
        if P(i,j) == d+1 
            C = kron(C, eye(n));
        else
            C = kron(C, A(:,P(i,j)));
        end
    end
    B = B + C;
end