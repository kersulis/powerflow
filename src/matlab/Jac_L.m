function L = Jac_L(i,k,theta,V,Y,Q_i)
% This function computes an upper-left Jacobian element "H" given angles,
% voltage, admittance, and reactive power.
G = real(Y);
B = imag(Y);

if i == k
    L = -B(i,i)*(V(i)^2) + Q_i;
else
    L = V(i)*V(k)*(G(i,k)*sin(theta(i)-theta(k)) - B(i,k)*cos(theta(i)-theta(k)));
end