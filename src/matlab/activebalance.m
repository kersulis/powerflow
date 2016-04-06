function P = activebalance(i,theta,V,Y,n)
% This function accepts a row index i and two n-by-1 vectors theta and V.
% It computes active power using the active power balance equation (number
% 1 in Power Flow Equations handout).
G = real(Y);
B = imag(Y);
P = 0;

for k = 1:n
    P = P + V(k)*(G(i,k)*cos(theta(i)-theta(k)) + ...
                      B(i,k)*sin(theta(i)-theta(k)));
end
P = P*V(i);
end
    