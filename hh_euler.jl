using Plots, BenchmarkTools, Distributions

function euler(f::Function, u0::Vector{Float64}, p::Vector{Float64},tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(100/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 1;
    for i in Int(100/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    return solution(t,u)
end

struct solution
    t::Vector{Float64}
    u::Matrix{Float64}
end

# Gate functions

αn(V)=(0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0))
αm(V)=(0.182*(V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0))
αh(V)= 0.25 * exp((-1.0 * (V + 90.0)) / 12.0)

βn(V)=(-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0))
βm(V)=(-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0))
βh(V)=(0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)


function hodg_hux_det(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables
    V = u[1]
    n = u[2]
    m = u[3]
    h = u[4]

    # Channel currents
    I_na =  g_na * m^3 * h * (V - V_na)
    I_k  =  g_k * n^4 * (V - V_k)
    I_l  =  g_l * (V- V_l)
   
    # ODE system
     dV =  1/C * (I_tot -I_na - I_k - I_l)
     dn =  αn(V) * (1 - n) - βn(V)*n
     dm =  αm(V) * (1 - m) - βm(V)*m
     dh =  αh(V) * (1 - h) - βh(V)*h
    return [dV,dn,dm,dh]
end

# Parameters

V_na = 55.0;
V_k  = -77.0; 
V_l  = -65.0;
g_na = 40.0;
g_k  = 35.0;
g_l  = 0.3;
C = 1.0;
I_tot = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot]

# Initial conditions
n_inf(v) = αn(v) / (αn(v) + βn(v))
m_inf(v) = αm(v) / (αm(v) + βm(v))
h_inf(v) = αh(v) / (αh(v) + βh(v))
v₀ = -60.0
u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)]

# Integration with euler's method

tspan = (0,1000);
h = 1e-3;
sol = euler(hodg_hux_det, u₀, p, tspan, h)

# Plots

plot(sol.t,sol.u[1,:])
plot(sol.t,sol.u[2,:])
plot!(sol.t,sol.u[3,:])
plot!(sol.t,sol.u[4,:])