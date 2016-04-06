function Q = reactivebalance(i,theta,V,Y,n)
% This function accepts a row index i and two n-by-1 vectors theta and V.
% It computes reactive power using the reactive power balance equation (2
% in Power Flow Equations handout).
G = real(Y);
B = imag(Y);

Q = 0;

for k = 1:n
    Q = Q + V(k)*(G(i,k)*sin(theta(i)-theta(k)) - ...
                  B(i,k)*cos(theta(i)-theta(k)));
end
Q = Q*V(i);
end