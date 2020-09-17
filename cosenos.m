function C = cosenos(A,phi,N)

C = zeros(2*N,length(A));
for i=1:N
    C(i,:) = A.*cos(i*phi);
end;

for i=N+1:2*N
    C(i,:) = A.*sin((i-N)*phi);
end;