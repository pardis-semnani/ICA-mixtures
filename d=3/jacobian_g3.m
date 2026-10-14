% Finds the Jacobian of Q_{n,1,3} (1 component and d=3) at a give parameter
% (b, A, D, b1), where b: n*1, A: n*n, D: 1*n, b1: n*1, and b1 is the b 
% in the first component.
function J = jacobian_g3(b, A, D, b1) 

n = size(b,1);
J = zeros(n+n^2+n^3, 2*n+n^2+1);
C = 1/2 * der([b b]);
for j=1:n
    C = C + 1/2 * der([A(:,j) A(:,j)]);
end
J(:,1:n) = [eye(n) ; der(b) ; C];
for j=1:n
    J(:,j*n+1:(j+1)*n) = [zeros(n,n) ; der(A(:,j)) ; 
        1/2 * D(1,j) * der([A(:,j) A(:,j)]) + der([b A(:,j)])];
end
for j=1:n
    J(:,n^2+n+j) = [zeros(n,1) ; zeros(n^2,1) ; 
        kron(kron(A(:,j),A(:,j)),A(:,j))];
end
J(:, end) = [b - b1; kron(b,b) - kron(b1,b1) ; kron(kron(b,b),b) - kron(kron(b1,b1),b1)];
end