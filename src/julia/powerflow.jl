"""
Hiskens power flow translated to Julia.

IN

* `V` [pu], Vector of voltage magnitudes
* `T` [rad], Vector of voltage angles
* `P` [pu], Vector of active power injections
* `Q` [pu], Vector of reactive power injections
* `ty`[-], Vector of bus types (1: PQ, 2: PV, 3: slack)
* `Vmax` [pu], Vector of maximum voltage magnitudes for PQ buses
* `Y` [pu], Admittance matrix
* `tol` [-], Newton iteration tolerance
* `maxiter` [-], Maximum allowable iterations

OUT

* `V` [pu], Vector of converged voltage magnitudes
* `T` [pu], Vector of converged voltage angles
* `Pcalc` [pu], Vector of calculated active power injections
* `Qcalc` [pu], Vector of calculated reactive power injections
* `ty` [-], Vector of bus types

_Note: if a PQ bus reaches its maximum voltage, the algorithm converts it
to a PV bus. Thus, output `ty` may differ from input `ty`._
"""
function powerflow(
    V::Vector{Float64},
    T::Vector{Float64},
    P::Vector{Float64},
    Q::Vector{Float64},
    ty::Vector{Int64},
    Vmax::Vector{Float64},
    Y::AbstractArray;
    tol=1e-5,
    maxiter=20
    )

    G, B = real(Y), imag(Y)
    nbus = length(ty)

    converged = false
    Pcalc = zeros(nbus)
    Qcalc = zeros(nbus)

    for iter = 1:maxiter
        Pcalc = zeros(nbus)
        Qcalc = zeros(nbus)

        for i = 1:nbus
            for k in find(Y[i,:])
                Tik = T[i] - T[k]
                Pcalc[i] += V[i]*V[k]*(G[i,k]*cos(Tik) + B[i,k]*sin(Tik))
                Qcalc[i] += V[i]*V[k]*(G[i,k]*sin(Tik) - B[i,k]*cos(Tik))
            end

            if ty[i] == 2 && Vmax[i] > 0
                if Qcalc[i] > Q[i]
                    ty[i] = 1
                end
            end
        end

        dPda = spzeros(nbus, nbus)
        dPdV = spzeros(nbus, nbus)
        dQda = spzeros(nbus, nbus)
        dQdV = spzeros(nbus, nbus)
        Pmis = zeros(nbus)
        Qmis = zeros(nbus)

        for i = 1:nbus
            Vi = V[i]
            for k = find(Y[i,:])
                if i == k
                    dPda[i,i] = -B[i,i]*Vi*Vi - Qcalc[i]
                    dPdV[i,i] = G[i,i]*Vi + Pcalc[i]/Vi
                    dQda[i,i] = -G[i,i]*Vi*Vi + Pcalc[i]
                    dQdV[i,i] = -B[i,i]*Vi + Qcalc[i]/Vi
                else
                    Vk = V[k]
                    Tik = T[i] - T[k]
                    dPda[i,k] = Vi*Vk*(G[i,k]*sin(Tik) - B[i,k]*cos(Tik))
                    dPdV[i,k] = Vi*(G[i,k]*cos(Tik) + B[i,k]*sin(Tik))
                    dQda[i,k] = -Vi*Vk*(G[i,k]*cos(Tik) + B[i,k]*sin(Tik))
                    dQdV[i,k] = Vi*(G[i,k]*sin(Tik) - B[i,k]*cos(Tik))
                end
            end
        end

        for i=1:nbus
            if ty[i] == 1 # PQ bus
                Pmis[i] = Pcalc[i] - P[i]
                Qmis[i] = Qcalc[i] - Q[i]
            else # PV bus
                Pmis[i] = Pcalc[i] - P[i]
                Qmis[i] = 0
                dQda[i,:] = zeros(nbus)
                dQdV[i,:] = zeros(nbus)
                dPdV[:,i] = zeros(nbus) # PV bus: V constant
                dQdV[:,i] = zeros(nbus)
                dQdV[i,i] = 1
            end

            if ty[i] == 3 # slack bus
                Pmis[i] = 0
                dPda[i,:] = zeros(nbus)
                dPdV[i,:] = zeros(nbus)
                dPda[:,i] = zeros(nbus)
                dQda[:,i] = zeros(nbus)
                dPda[i,i] = 1
            end
        end

        mis = [Pmis;Qmis]

        if maxabs(mis)<tol
            converged = true
            info("Converged, $(iter - 1) iterations")
            break
        else
            J = [dPda dPdV;dQda dQdV]
            update = -J\mis #full operator and result (not sparse)

            T += update[1:nbus]
            V += update[nbus+1:2*nbus]

            # fix LTC-regulated voltages
            # for t in LTC
            #     V[t.reg] = V[t.obs] + t.alpha*t.tap
            # end

            for i=1:nbus
                if ty[i] == 1 && Vmax[i] > 0
                    if V[i] > Vmax[i]
                        V[i] = Vmax[i]
                        ty[i] = 2 # switch from PQ to PV bus
                        info("Bus $i has hit voltage limit when Q = $(Q[i])")
                    end
                end 
            end
        end
    end
    !converged && warn("Powerflow not converged")
    return V,T,Pcalc,Qcalc,ty
end

"""
    (Y,Yf,Yt) = getY(lineOut,f,t,r,x,b,tap,ysh)
This function builds admittance matrices.

INPUTS:

*  `lineOut`: [nline x 1, logical] which lines are out of service
*  `f`: nline Vector of from buses (range 1 to nbus)
*  `t`: nline Vector of to buses (range 1 to nbus)
*  `r`: nline Vector of series resistances (pu)
*  `x`: nline Vector of series reactances (pu)
*  `b`: nline Vector of line charging shunt susceptances (pu)
*  `tap`: nline vector of complex tap ratios (pu)
*  `ysh`: nbus vector of complex bus shunt admittances (pu)

for each branch, compute the elements of the branch admittance matrix where
```
  | If |   | Yff  Yft |   | Vf |
  |    | = |          | * |    |
  | It |   | Ytf  Ytt |   | Vt |
```
*Assumes the grid is connected (i.e. no island nodes).*
"""
function getY(
    from::Vector{Int64},
    to::Vector{Int64},
    r::Vector{Float64},
    x::Vector{Float64},
    b::Vector{Float64},
    tap::Vector{Complex{Float64}} = fill(Complex(1.), length(from)),
    ysh::Vector{Complex{Float64}} = fill(0.im, length(unique([from;to]))),
    lineOut::Vector{Bool} = fill(false, length(from))
    )
    # remap bus tags to 1:nbus indices
    tags = sort(unique([from;to]))
    f = tags2indices(from,tags)
    t = tags2indices(to,tags)

    inService = find(!lineOut)
    tap = tap[inService]
    fis = f[inService]
    tis = t[inService]

    Ys = (1./(r + x*im))[inService]
    Bc = b[inService]
    Ytt = Ys + Bc*im/2
    Yff = Ytt./(tap.*conj(tap))
    Yft = -Ys./conj(tap)
    Ytf = -Ys./tap

    # build connection matrices
    nbus,nline = length(ysh),length(lineOut)
    # connection matrix for line and from buses
    Cf = sparse(1:nline,f,ones(nline),nline,nbus)
    # connection matrix for line and to buses
    Ct = sparse(1:nline,t,ones(nline),nline,nbus)

    # build Yf and Yt such that Yf*V is the vector of complex branch currents injected.
    i = [inService;inService]
    Yf = sparse(i,[fis;tis],[Yff;Yft],nline,nbus)
    Yt = sparse(i,[fis;tis],[Ytf;Ytt],nline,nbus)

    # build Ybus
    Y = Cf'*Yf + Ct'*Yt + sparse(1:nbus,1:nbus,ysh,nbus,nbus)
    return Y #,Yf,Yt
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
