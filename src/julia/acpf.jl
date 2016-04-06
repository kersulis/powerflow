type Bus
    "bus type"
    ty::ASCIIString
    "Net active power [pu]"
    P::Float64
    "Net reactive power [pu]"
    Q::Float64
    "Voltage magnitude [pu]"
    V::Float64
    "Voltage angle [rad]"
    T::Float64
end

type Line
    "from"
    f::Int64
    "to"
    t::Int64
    "resistance"
    r::Float64
    "reactance"
    x::Float64
    "susceptance"
    b::Float64
end

"""
Also known as dPdT.
"""
function Hik(i, k, Y, V, T, P, Q)
    if i == k
        Bii = imag(Y[i,i])
        return -Bii*(V[i]^2) - Q[i]
    else
        Gik, Bik = real(Y[i,k]), imag(Y[i,k])
        Tik = T[i] - T[k]
        return V[i]*V[k]*(Gik*sin(Tik) - Bik*cos(Tik))
    end
end

"""
Also known as dPdV.
"""
function Nik(i, k, Y, V, T, P, Q)
    if i == k
        Gii = real(Y[i,i])
        return Gii*V[i] + P[i]/V[i]
    else
        Gik, Bik = real(Y[i,k]), imag(Y[i,k])
        Tik = T[i] - T[k]
        return V[i]*(Gik*cos(Tik) + Bik*sin(Tik))
    end
end

"""
Also known as dQdT.
"""
function Jik(i, k, Y, V, T, P, Q)
    if i == k
        Gii = real(Y[i,i])
        return -Gii*(V[i]^2) + P[i]
    else
        return -Nik(i, k, Y, V, T, P, Q)
    end
end

"""
Also known as L_ik.
"""
function Lik(i, k, Y, V, T, P, Q)
    if i == k
        Bii = imag(Y[i,i])
        return -Bii*V[i] + Q[i]/V[i]
    else
        Gik, Bik = real(Y[i,k]), imag(Y[i,k])
        Tik = T[i] - T[k]
        return V[i]*(Gik*sin(Tik) - Bik*cos(Tik))
    end
end

"""
This function computes all rows and columns of the Jacobian. When using
this function in a power flow algorith, be sure to use only those rows
and columns that correspond to update (unknown) quantities.
"""
function pf_jacobian(buses, Y, V, T, P, Q)
    p_idx = find([b.ty != "slack" for b in buses])
    q_idx = find([b.ty == "PQ" for b in buses])

    H = [Float64(Hik(i, k, Y, V, T, P, Q)) for i in p_idx, k in p_idx]
    N = [Float64(Nik(i, k, Y, V, T, P, Q)) for i in p_idx, k in q_idx]
    J = [Float64(Jik(i, k, Y, V, T, P, Q)) for i in q_idx, k in p_idx]
    L = [Float64(Lik(i, k, Y, V, T, P, Q)) for i in q_idx, k in q_idx]
    return [H N;J L]
end

"""
This function accepts a row index i and two n-by-1 vectors theta and V.
It computes active mismatch using the active power balance equation (number
1 in Power Flow Equations handout).
"""
function mismatch_p(buses, Y, V, T)
    n = length(buses)
    p_idx = 1:n #find([b.ty != "slack" for b in buses])
    Pvec = Vector{Float64}()
    for i in p_idx
        P = 0.0
        for k = 1:n
            Gik, Bik = real(Y[i,k]), imag(Y[i,k])
            Tik = T[i] - T[k]
            P += V[k]*(Gik*cos(Tik) + Bik*sin(Tik))
        end
        push!(Pvec, buses[i].P - P*V[i])
    end
    return Pvec
end

"""
This function accepts a row index i and two n-by-1 vectors theta and V.
It computes reactive mismatch using the reactive power balance equation (2
in Power Flow Equations handout).
"""
function mismatch_q(buses, Y, V, T)
    n = length(buses)
    q_idx = 1:n #find([b.ty == "PQ" for b in buses])
    Qvec = Vector{Float64}()
    for i in q_idx
        Q = 0.0
        for k = 1:n
            Gik, Bik = real(Y[i,k]), imag(Y[i,k])
            Tik = T[i] - T[k]
            Q += V[k]*(Gik*sin(Tik) - Bik*cos(Tik))
        end
        push!(Qvec, buses[i].Q - Q*V[i])
    end
    return Qvec
end

"""
    createY(f,t,x [,r,b]) -> Y
Create an admittance matrix for AC power flow.
All inputs are real. The output matrix is real if no line
resistances are provided (DC case), and complex otherwise.
* `f`,`t`: vectors encoding all lines (fi,ti)
* `x`: per-unit reactance xi for all lines
* `r`: per-unit resistance ri for all lines
* `b`: per-unit susceptance bi for all lines
"""
function createY(
    f::Vector{Int64},
    t::Vector{Int64},
    x::Vector{Float64},
    r=0.0::Union{Vector{Float64},Float64},
    b=0.0::Union{Vector{Float64},Float64}
    )
    z = r + x*1im
    y = 1./z
    b = b*1im
    Y = sparse([f; t; t; f],[t; f; t; f],[-y; -y; y + b./2; y + b./2])

    # for DC power flow, we typically want a matrix with real entries:
    if r == 0
        return imag(Y)
    else
        return Y
    end
end

function createY(lines::Vector{Line})
    f = [Int64(l.f) for l in lines]
    t = [Int64(l.t) for l in lines]
    r = [Float64(l.r) for l in lines]
    x = [Float64(l.x) for l in lines]
    b = [Float64(l.b) for l in lines]
    return createY(f, t, x, r, b)
end

function acpf(
    buses::Vector{Bus},
    lines,
    max_iter,
    tol
    )
    n = length(buses)
    # slack_idx = findfirst([b.ty == "slack" for b in buses])

    Y = createY(lines)
    V = [b.V for b in buses]
    T = [b.T for b in buses]
    p_idx = find([b.ty != "slack" for b in buses])
    q_idx = find([b.ty == "PQ" for b in buses])

    iter = 0
    converged = false
    P = zeros(n)
    Q = zeros(n)
    while iter <= max_iter && !converged
        P = mismatch_p(buses, Y, V, T)
        Q = mismatch_q(buses, Y, V, T)
        converged = maxabs([P[p_idx];Q[q_idx]]) < tol
        # println(iter)
        if converged
            info("Converged after $iter iterations")
        end
        iter += 1
        J = pf_jacobian(buses, Y, V, T, P, Q)
        update = J\[P[p_idx];Q[q_idx]]
        T[p_idx] += update[1:length(p_idx)]
        V[q_idx] += update[(end-length(q_idx)+1):end]
    end
    return P, Q, V, T
end
