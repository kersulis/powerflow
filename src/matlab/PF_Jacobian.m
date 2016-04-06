function J = PF_Jacobian(n,Y,P,Q,theta,V)
% This function computes all rows and columns of the Jacobian. When using
% this function in a power flow algorith, be sure to use only those rows
% and columns that correspond to update (unknown) quantities.

G = real(Y);
B = imag(Y);
J = zeros(n);

% H, the upper-left block:
for r = 1:n
    for c = 1:n
        if r == c
            J(r,c) = -B(r,c)*(V(r)^2) - Q(r);
        else
            J(r,c) = V(r)*V(c)*(G(r,c)*sin(theta(r)-theta(c)) - ...
                     B(r,c)*cos(theta(r)-theta(c)));
        end
    end
end

% N, the upper-right block:
for r = 1:n
    for c = (n + 1):(n*2)
        if r == (c-n)
            J(r,c) = G(r,c-n)*(V(r)^2) + P(r);
        else
            J(r,c) = V(r)*V(c-n)*(G(r,c-n)*cos(theta(r)-theta(c-n)) + ...
                     B(r,c-n)*sin(theta(r)-theta(c-n)));
        end
    end
end

% J, the lower-left block:
for r = (n + 1):(n*2)
    for c = 1:n
        if (r-n) == c
            J(r,c) = -G(r-n,c)*(V(r-n)^2) + P(r-n);
        else
            J(r,c) = -V(r-n)*V(c)*(G(r-n,c)*cos(theta(r-n)-theta(c)) + ...
                     B(r-n,c)*sin(theta(r-n)-theta(c)));
        end
    end
end

% L, the lower-right block:
for r = (n + 1):(n*2)
    for c = (n + 1):(n*2)
        if r == c
            J(r,c) = -B(r-n,c-n)*(V(r-n)^2) + Q(r-n);
        else
          J(r,c)=V(r-n)*V(c-n)*(G(r-n,c-n)*sin(theta(r-n)-theta(c-n)) - ...
                     B(r-n,c-n)*cos(theta(r-n)-theta(c-n)));
        end
    end
end

end
