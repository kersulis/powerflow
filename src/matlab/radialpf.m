function [P, Q, V] = radialpf(nodes,LD,pL,qL,LTCnodes,LTCtaps)
% Radial Power Flow
% Jonas Kersulis
% This code runs a power flow on a radial system.
% 'nodes' is a list of unique nodes sorted so that every node below a
% certain index is "downstream".
% 'network' is a matrix with a 1 wherever two nodes are connected.
% 'LD' is a struct with line data. It has fields for NodeA (from), NodeB
% (to), R (resistance), and X (reactance).
% 'pL' and 'qL' are active and reactive load at each node, sorted by 'nodes'.
% 'LTCnodes' is a vector where each element is a node that is downstream
% from a tap-changing transformer. 'LTCtaps' is a vector with the
% corresponding turns ratio for each transformer.

r = LD.R;
x = LD.X;
NodeA = LD.NodeA;
NodeB = LD.NodeB;

% Switch to internal numbering:
for i = 1:length(NodeA)
    for j = 1:length(nodes)
        if nodes(j) == NodeA(i)
            NodeA(i) = j;
        end
        if nodes(j) == NodeB(i)
            NodeB(i) = j;
        end
    end
end

% 1. Initialize all power flows to 0
P = zeros(length(nodes),1);
Q = P;
sending = zeros(32,1);

% 2. Initialize all voltages to 1.05
V(1:length(nodes),1) = 1.05;

% Find dead ends:
for i = 2:length(nodes)
    % Find sending node:
    temp = find(NodeA == i);
    adj = [NodeB(temp) r(temp) x(temp)];
    temp = find(NodeB == i);
    adj = [adj; [NodeA(temp) r(temp) x(temp)]];   % adj is connected nodes.
    adj = adj(adj(:,1) < i,:);  % adj is sending node for node i.
    sending(i) = adj(1);% Sending node for node i ("i-1")
    r(i-1) = adj(2);    % From sending to receiving (i-1 to i)
    x(i-1) = adj(3);    % From sending to receiving (i-1 to i)
end
% Nodes that do not send to other nodes are ends. Should be nodes 5, 12,
% 14, 18, 22, 27, 30, 32.
endsmaster = setdiff(1:length(nodes),sending);
counter = 0;
converged = 0;

% Iterate power and voltage sweeps until convergence is reached.
while converged == 0
    Ptest = P; % Used later to test for convergence
    ends = endsmaster; % Refresh ends (empty vector after each sweep)

    % Injections at tips are just loads:
    for i = 1:length(ends)
        P(ends(i)) = pL(ends(i));
        Q(ends(i)) = qL(ends(i));
    end

    % POWER SWEEP: Start at tip furthest from feeder (node 32).
    % Use load at tip to estimate injection into next node.
    tip = ends(end);
    j = tip; % Initialize j to node furthest from feeder (32)
    ends = setdiff(ends,tip); % remove current tip from ends
    tip = ends(end);

    fork = 0; Pfork = 0; Qfork = 0; % Initialize some variables

    while j > 1 % Power sweep ends when feeder node injections are calculated.
        % Want: P(i)
        % Have: P(j)
        i = sending(j); % i is the node feeding node j

        % Store old injections into Pi and Qi:
        Pi = P(i); Qi = Q(i);

        % Store recently calculated injections in Pj and Qj:
        Pj = P(j); Qj = Q(j);

        % Store voltage at node i:
        Vi = V(i);

        % Store resistance and reactance between i and j:
        rij = r(j-1); xij = x(j-1);

        % Store real and reactive load at j:
        pLj = pL(j); qLj = qL(j);

        % Calculate active and reactive injections at node i:
        P(i) = Pj + rij*(Pi^2 + Qi^2)/(Vi^2) + pLj;
        Q(i) = Qj + xij*(Pi^2 + Qi^2)/(Vi^2) + qLj;

        % This completes power injection calculations for a typical node. But
        % what if there are two branches leaving node i?
        if i < tip % True if there is a second branch to jump ahead to
            j = tip; % Jump ahead to tip
            ends = setdiff(ends,tip); % Remove new tip from ends

            if isempty(ends) % All out of ends
                tip = 0;
            else
                tip = ends(end); % Store new tip
            end

            % Temporarily store values at junction to add later:
            fork = i;
            Pfork = P(i);
            Qfork = Q(i);
        else
            j = i; % j becomes i
        end

        % If we are reaching junction node i from a second branch, we need to
        % add injections calculated from first branch (Pfork and Qfork):
        if i == fork
            P(i) = P(i) + Pfork;
            Q(i) = Q(i) + Qfork;
        end
    end

    % VOLTAGE SWEEP: Start at feeder node
    for j = 2:length(nodes)
        i = sending(j);
        % Want: Vj
        % Have: Vi
       Vi = V(i);

       rij = r(j-1); xij = x(j-1);

       Pi = P(i); Qi = Q(i);

       Vj = Vi^2 - 2*(rij*Pi + xij*Qi) + (rij^2 + xij^2)*((Pi^2 + Qi^2)/(Vi^2));
       V(j) = sqrt(Vj);

       % Fix voltages downstream from tap-changing transformers:
       index = find(LTCnodes == j);
       if any(index)
           V(j) = LTCtaps(index);
       end
    end

    % At this point, a complete power sweep/voltage sweep has been performed.
    % To determine convergence, take the difference between updated values and
    % past values.

    counter = counter + 1;
    if abs(max(P - Ptest)) < 1e-5
        converged = 1;
    end
end
end
