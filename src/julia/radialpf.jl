"""
Radial Power Flow
Jonas Kersulis, 2016-03-21

Run [DistFlow](http://dx.doi.org/10.1109/61.19266) on a radial power system.

IN:

* 'nodes', a dictionary with keys `index` (unique nodes sorted such that 
every node below a certain node's index is downstream of that node), `pL` (active load at each node), and `qL` (reactive load at each node).
* 'network', a matrix with a 1 wherever two nodes are connected.
* 'lines', a line data dictionary with keys `f` (from), `t`
(to), `r` (resistance), and `x` (reactance).
* 'LTCnodes', a vector where each element is a node that is downstream
from a tap-changing transformer. 
* 'LTCtaps', a vector with the corresponding turns ratio for each transformer.

OUT:

* `P`: active power flows across all lines.
* `Q`: reactive power flows across all lines.
* `V`: voltage magnitudes at all nodes.
"""
function radialpf(nodes,lines,LTCnodes,LTCtaps)
    index = nodes["index"]
    pL = nodes["pL"]
    qL = nodes["qL"]
    nb = length(index) # number of buses

    r = lines["r"]
    x = lines["x"]
    f = lines["f"]
    t = lines["t"]


    # Switch to internal numbering:
    f = tags2indices(f, index)
    t = tags2indices(t, index)

    # 1. Initialize all power flows to 0
    P = zeros(nb)
    Q = zeros(nb)
    sending = zeros(nb)

    # 2. Initialize all voltages to 1.0
    V[1:nb,1] = 1.0

    # Find dead ends:
    sending = [0]
    for i = 2:nb
        # Find sending node:
        temp = find(f .== i)
        adj = [t[temp] r[temp] x[temp]]
        temp = find(t .== i)
        adj = [adj; [f[temp] r[temp] x[temp]]] # adj is connected nodes.
        adj = adj[adj[:,1] .< i,:]             # adj is sending node for node i.
        push!(sending, adj[1])                 # Sending node for node i ("i-1")
        r[i-1] = adj[2]                        # From sending to receiving (to i)
        x[i-1] = adj[3]                        # From sending to receiving (to i)
    end
    # Nodes that do not send to other nodes are ends. 
    endsmaster = setdiff(1:nb,sending)
    counter = 0
    converged = false

    # Iterate power and voltage sweeps until convergence is reached.
    while !converged
        Ptest = P         # Used later to test for convergence
        ends = endsmaster # Refresh ends (empty vector after each sweep)

        # Injections at tips are just loads:
        for i in ends
            P[i] = pL[i]
            Q[i] = qL[i]
        end
        # POWER SWEEP: Start at tip furthest from feeder (node 32).
        # Use load at tip to estimate injection into next node.
        tip = ends[end]
        j = tip                  # Initialize j to node furthest from feeder (32)
        ends = setdiff(ends,tip) # remove current tip from ends
        tip = ends[end]

        # Initialize some variables
        fork = 0
        Pfork = 0.0
        Qfork = 0.0

        while j > 1 # Power sweep ends when feeder node injections are calculated.
            # Want: P(i)
            # Have: P(j)
            i = sending[j] # i is the node feeding node j

            # Store old injections into Pi and Qi:
            Pi = P[i]
            Qi = Q[i]

            # Store recently calculated injections in Pj and Qj:
            Pj = P[j]
            Qj = Q[j]

            # Store voltage at node i:
            Vi = V[i]

            # Store resistance and reactance between i and j:
            rij = r[j-1]
            xij = x[j-1]

            # Store real and reactive load at j:
            pLj = pL[j]
            qLj = qL[j]

            # Calculate active and reactive injections at node i:
            P[i] = Pj + rij*(Pi^2 + Qi^2)/(Vi^2) + pLj
            Q[i] = Qj + xij*(Pi^2 + Qi^2)/(Vi^2) + qLj

            # This completes power injection calculations for a typical node. But
            # what if there are two branches leaving node i?
            if i < tip  # True if there is a second branch to jump ahead to
                j = tip # Jump ahead to tip
                ends = setdiff(ends,tip) # Remove new tip from ends

                if isempty(ends) # All out of ends
                    tip = 0
                else
                    tip = ends[end] # Store new tip
                end

                # Temporarily store values at junction to add later:
                fork = i
                Pfork = P[i]
                Qfork = Q[i]
            else
                j = i # j becomes i
            end

            # If we are reaching junction node i from a second branch, we need to
            # add injections calculated from first branch (Pfork and Qfork):
            if i == fork
                P[i] = P[i] + Pfork
                Q[i] = Q[i] + Qfork
            end
        end

        # VOLTAGE SWEEP: Start at feeder node
        for j = 2:nb
            i = sending[j]
            # Want: Vj
            # Have: Vi
           Vi = V[i]

           rij = r[j-1]
           xij = x[j-1]

           Pi = P[i]
           Qi = Q[i]

           Vj = Vi^2 - 2*(rij*Pi + xij*Qi) + (rij^2 + xij^2)*((Pi^2 + Qi^2)/(Vi^2))
           V[j] = sqrt(Vj)

           # Fix voltages downstream from tap-changing transformers:
           LTCidx = find(LTCnodes .== j)
           if !isempty(LTCidx)
               V[j] = LTCtaps[LTCidx]
           end
        end

        # At this point, a complete power sweep/voltage sweep has been performed.
        # To determine convergence, take the difference between updated values and
        # past values.

        counter += 1
        if maxabs(P - Ptest) < 1e-5
            converged = true
        end
    end
end

"""
Given a vector `v` containing bus tags, and a complete and sorted
vector of bus tags, return a version of `v` mapped to bus indices.

Example: if there are three buses in the network called 101, 202, and 303,
one would let `tags = [101;202;303]`. Corresponding indices are simply `[1;2;3]`. So this function would convert a vector `v = [202;303]` into `vidx = [2;3]`.
"""
function tags2indices(
    v::Vector{Int64},
    tags::Vector{Int64}
    )
    indices = collect(1:length(tags))
    vidx = Vector{Int64}()
    for vi in v
        push!(vidx,indices[findfirst(tags,vi)])
    end
    return vidx
end

