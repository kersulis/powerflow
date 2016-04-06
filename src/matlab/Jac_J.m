function J = Jac_J(i,k,theta,V,Y,P_i)
% This function computes an upper-right Jacobian element "H" given angles,
% voltage, admittance, and active power.
G = real(Y);
B = imag(Y);

if i == k
    J = -G(i,i)*(V(i)^2) + P_i;
else
    J = -V(i)*V(k)*(G(i,k)*cos(theta(i)-theta(k)) + B(i,k)*sin(theta(i)-theta(k)));
end