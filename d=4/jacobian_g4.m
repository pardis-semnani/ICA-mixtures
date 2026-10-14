% Finds the Jacobian of Q_{n,1,4} (1 component and d=4) 
% at a give parameter (b, A, D, F, b1), where b: n*1, A: n*n, D: 1*n, F: 1*n, 
% b1: n*1, and b1 is the b in the first component.
function J = jacobian_g4(b, A, D, F, b1) 
n = size(b,1);
J = zeros(n+n^2+n^3+n^4, 3*n+n^2+1);
C = 1/2 * der([b b]);
C1 = 1/6 * der([b b b]);
for j=1:n
    C = C + 1/2 * der([A(:,j) A(:,j)]);
    C1 = C1 + 1/6 * D(1,j) * der([A(:,j) A(:,j) A(:,j)]) + 1/2 * der([b A(:,j) A(:,j)]);
end
J(:,1:n) = [eye(n) ; der(b) ; C ; C1];
for j=1:n
    C2 = (F(1,j) + 1.5)/6 * der([A(:,j) A(:,j) A(:,j)]) + D(1,j)/2 * der([b A(:,j) A(:,j)]) + 1/2 * der([A(:,j) b b]);
    for l = 1:n
        C2 = C2 + 1/4 * der([A(:,j) A(:,l) A(:,l)]);
    end
    J(:,j*n+1:(j+1)*n) = [zeros(n,n) ; der(A(:,j)) ; 
        1/2 * D(1,j) * der([A(:,j) A(:,j)]) + der([b A(:,j)]) ; C2];
end
for j=1:n
    J(:,n^2+n+j) = [zeros(n,1) ; zeros(n^2,1) ; 
        kron(kron(A(:,j),A(:,j)),A(:,j)) ;
        kron(kron(kron(A(:,j),A(:,j)),A(:,j)),b) + kron(kron(kron(A(:,j),A(:,j)),b),A(:,j)) + kron(kron(kron(A(:,j),b),A(:,j)),A(:,j)) + kron(kron(kron(b,A(:,j)),A(:,j)),A(:,j))];
end
for j = 1:n
    J(:, n^2+2*n+j) = [zeros(n+n^2+n^3,1) ; kron(kron(kron(A(:,j),A(:,j)),A(:,j)),A(:,j))];
end

J(:, end) = [b - b1; kron(b,b) - kron(b1,b1) ; kron(kron(b,b),b) - kron(kron(b1,b1),b1) ;
    kron(kron(kron(b,b),b),b) - kron(kron(kron(b1,b1),b1),b1)];
end