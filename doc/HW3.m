%%   EECS 598 Homework 3
%   Jonas Kersulis

% Part 1
% In Part 1, bus 4 is a PQ bus--its active and reactive power are specified
% to be zero, and both mismatches must be calculated in each iteration to
% solve for angle and voltage.

%% Step 1: Find Admittance Matrix
clc
clear all
% Translate network topology and parameters to matrix form by
% building the admittance matrix Y.

%   Bus types:
%   0   Slack
%   1   PQ
%   2   PV

%        Num  Type  P       Q     theta      V  
BusInfo =  [1   0   0       0       0       1.02;
            2   2   0.4     0       0       1.05;
            3   1   -1.0    -0.3    0       0;
            4   1   0       0       0       0;];
        
n = size(BusInfo,1);               

%         From  To   R       X        B
LineInfo = [1   2   0.10    1.00    0.50;
            1   3   0.05    0.60    0.25;
            2   4   0.05    0.40    0.25;
            3   4   0.00    0.0001  0.00;];

% Mutual admittance is found by taking the inverse of the complex impedance
% between each connected pair of buses:
y_mutual = [LineInfo(:,1), LineInfo(:,2), (1./complex(LineInfo(:,3), LineInfo(:,4)))];

% Self admittance is the sum of the admittances of all lines connected to
% each bus, plus half of the total charging susceptance associated with
% those lines.
y_self =   [BusInfo(1,1), BusInfo(1,1), (y_mutual(1,3)+y_mutual(2,3)+ 0.5i.*LineInfo(1,5) + 0.5i.*LineInfo(2,5));
            BusInfo(2,1), BusInfo(2,1), (y_mutual(1,3)+y_mutual(3,3)+ 0.5i.*LineInfo(1,5) + 0.5i.*LineInfo(3,5));
            BusInfo(3,1), BusInfo(3,1), (y_mutual(2,3)+y_mutual(4,3)+ 0.5i.*LineInfo(2,5) + 0.5i.*LineInfo(4,5));
            BusInfo(4,1), BusInfo(4,1), (y_mutual(3,3)+y_mutual(4,3)+ 0.5i.*LineInfo(3,5) + 0.5i.*LineInfo(4,5));];

% The admittance matrix is found by vertically concatenating the mutual and
% self admittance matrices as follows:
Y = sparse([y_mutual(:,1); y_mutual(:,2); y_self(:,1)],...
           [y_mutual(:,2); y_mutual(:,1); y_self(:,2)],...
           [-1.*y_mutual(:,3); -1.*y_mutual(:,3); y_self(:,3)]);

%% Step 2: Initialize with Givens and Guesses

% Bus 1: slack bus (no equations)
% theta1    := 0
% V1        := 1.02pu

% Bus 2: PV bus
% Find active power mismatch
% P2_sp     := 0.4pu
% V2        := 1.05pu

% Bus 3: PQ bus
% Find active power mismatch
% Find reactive power mismatch
% be calculated)
% P3_sp     := -1.0pu
% Q3_sp     := -0.3pu

% Bus 4: PQ bus
% Find active power mismatch
% Find reactive power mismatch
% P4_sp     := 0
% Q4_sp     := 0

num_unknowns = 5;
% We need to make an initial guess for each of the five unknowns:
theta2_guess = 0.1;%input('Guess for theta2: ');
theta3_guess = 0.1;%input('Guess for theta3: ');
theta4_guess = 0.1;%input('Guess for theta4: ');
V3_guess = input('Guess for V3: ');
V4_guess = input('Guess for V4: ');

% Use guesses and given information to construct a vector for each variable
% type:
P_sp = [0; 0.4; -1.0; 0];
Q_sp = [0; 0; -0.3; 0];
theta = [0; theta2_guess; theta3_guess; theta4_guess];
V = [1.02; 1.05; V3_guess; V4_guess];

% Instantiate nett active and reactive power vectors:
P_nett(1:n) = 0;
Q_nett(1:n) = 0;

Jacobian = zeros(num_unknowns); % Fix the size of the Jacobian
mismatch(1:num_unknowns) = 0;
mismatch = mismatch';

%% Newton-Raphson Implmentation
tolerance = 1e-5;   % When updates are smaller than this value, stop
                    % iterating.
max_iterations = 10;% Do not perform any more iterations than this.

% Set up matrix to store results:
Result(1:num_unknowns+1,1:max_iterations) = 0;

for j = 1:max_iterations
    % Calculate nett active and reactive power at each bus:
    for x = 1:n
        P_nett(x) = activebalance(x,theta,V,Y,n);
        Q_nett(x) = reactivebalance(x,theta,V,Y,n);
    end
        
    % Calculate mismatches:
    mismatch(1) = P_sp(2) - P_nett(2);   % Bus 2 active
    mismatch(2) = P_sp(3) - P_nett(3);   % Bus 3 active
    mismatch(3) = P_sp(4) - P_nett(4);   % Bus 4 active
    
    mismatch(4) = Q_sp(3) - Q_nett(3); % Bus 3 reactive
    mismatch(5) = Q_sp(4) - Q_nett(4); % Bus 4 reactive
        
    % Use H and N functions to find four parts of the Jacobian matrix:
    H = Jacobian(1:3,1:3);
    for x = 2:4
        for y = 2:4
            H(x-1,y-1) = Jac_H(x,y,theta,V,Y,Q_nett(y));
        end
    end
    
    N = Jacobian(1:3,4:5);
    for x = 2:4
        for y = 3:4
            N(x-1,y-2) = Jac_N(x,y,theta,V,Y,P_nett(y));
        end
    end
    
    J = Jacobian(4:5,1:3);
    for x = 3:4
        for y = 2:4
            J(x-2,y-1) = -Jac_N(x,y,theta,V,Y,P_nett(y));
        end
    end
    
    L = Jacobian(4:5,4:5);
    for x = 3:4
        for y = 3:4
            L(x-2,y-2) = Jac_H(x,y,theta,V,Y,Q_nett(y));
        end
    end
    
    % Piece together Jacobian using four matrices calculated above:
    Jacobian = [H N; J L];

    % Solve for unknowns:
    update = Jacobian\mismatch;
    
    % Use updated values in next iteration:
    theta(2) = theta(2) + update(1);
    theta(3) = theta(3) + update(2);
    theta(4) = theta(4) + update(3);
    
    V(3) = V(3)*(1 + update(4));   % Multiply update by V3 because
                                   % algorithm solved for del(V3)/V3
    V(4) = V(4)*(1 + update(5));
    
    % Store results for later review:
    Result(1,j) = j;
    Result(2:num_unknowns+1,j) = update;
    if(all(abs(update)<tolerance))
        break
    end
end
disp(Result)
            
    