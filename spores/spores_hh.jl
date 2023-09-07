using Plots, Distributions, LaTeXStrings, DifferentialEquations,StaticArrays

# Parameters
const g_k = 30; # potassium channel conductance (min^-1)
const g_kq = 1; # potassium channel conductance with quinine (min^-1)
const g_n = 0.08; # nonspecific channel conductance (min^-1)
const g_nq = 0.005; # nonspecific channel conductance with quinine (min^-1)
const V_k0 = 30; # Nernst potential prefactor (mV)
const V_n = -130; # nonspecific ion Nernst potential (mV)
const alpha_g = 3; # channel opening constant due to germinant (min^-1)
const beta = 0.6; # channel opening rate decay constant (min^-1)
const V_0wt = -79; # initial membrane potential, wild-type (mV)
const V_0ktrc = -81; # initial membrane potential, ΔktrC (mV)
const V_0yugO = -80; # initial membrane potential, ΔyugO (mV)
const gamma_e = 1; # extracellular potassium relaxation rate (min^-1)
const F = 5.6; # membrane capacitance (mM/mV)
const K_m = 8; # external media potassium (mM)
const K_s = 235; # potassium threshold triggering germination (mM)
const K_wt = 300; # average initial potassium concentration, wild-type (mM)
const K_ktrc = 275; # average initial potassium concentration, ΔktrC (mM)
const K_yugO = 280; # average initial potassium concentration, ΔyugO (mM)
const sigma_wt = 15; # st. dev. of initial potassium concentration, wild type (mM)
const sigma_ktrc = 15; # st. dev. of initial potassium concentration, ΔktrC (mM)
const sigma_yugO = 35; # st. dev. of initial potassium concentration, ΔyugO  (mM)
alpha = 0; # germinant not present
# alpha = alpha_g; # germinant present

p = [g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha];

function euler(f::Function, u0::Vector{Float64}, p::Vector{Float64},
    tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    p[22]=0;
    for i in 1:Int(3600/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=p[7];
    for i in Int(3601/h):Int(3781/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=0;
    for i in Int(3781/h):Int(10800/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=p[7];
    for i in Int(10801/h):Int(10980/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=0;
    for i in Int(10981/h):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end

    return solution_euler(t,u)
end

struct solution_euler
    t::Vector{Float64}
    u::Matrix{Float64}
end

# struct solution_bin
#     t::Vector{Float64}
#     V::Vector{Float64}
#     N4::Vector{Float64}
#     M3::Vector{Float64}
#     H::Vector{Float64}
#     intensitat::Vector{Float64}
# end

# struct solution_mar
#     t::Vector{Float64}
#     V::Vector{Float64}
#     N4::Vector{Float64}
#     M3::Vector{Float64}
#     H::Vector{Float64}

#     N0::Vector{Float64}
#     N1::Vector{Float64}
#     N2::Vector{Float64}
#     N3::Vector{Float64}
#     changes::Vector{Float64}
#     intensitat_vars::Vector{Float64}
# end
#-------------------------------------------------------------- deterministic, model
function spores_hh_det(u,p,t)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;
    V = u[1]
    K_e = u[2]
    K_i = u[3]
    n = u[4]

    # Nernst potential
    V_k = V_k0*log(abs(K_e/K_i))

    #Dynamical system
    dV = -g_k * n^4 * (V - V_k) - g_n * n^4 * (V - V_n)
    dK_e = F * g_k * n^4 * (V - V_k) + F * g_n * n^4 * (V - V_n) - gamma_e * (K_e - K_m)
    dK_i = -F * g_k * n^4 * (V - V_n) - F * g_n * n^4 * (V - V_n) 
    dn = alpha * (1 - n) - beta * n

    return [dV,dK_e,dK_i,dn]
end

# 15 hores = 54000 s
tspan = (0,14400);
h = 1e-3;
u₀_det=rand(4);
sol_det = euler(spores_hh_det, u₀_det, p, tspan, h);

# Plots
plot(sol_det.t,sol_det.u[1,:], label = "V")
plot(sol_det.t,sol_det.u[2,:],label="K_e")
plot!(sol_det.t,sol_det.u[3,:],label="K_i")
plot(sol_det.t,sol_det.u[4,:],label="n",ylabel="fraction of open subunits")

#---------------------------------------------------------------- deterministic, states
function hodg_hux_det_states(u, p, t)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO = p;
    V = u[1]
    K_e = u[2]
    K_i = u[3]
    n₀ = u[4]
    n₁ = u[5]
    n₂ = u[6]
    n₃ = u[7]
    n₄ = u[8]

    V_k = V_k0*log(abs(K_e/K_i))

    #Dynamical system
    dV = -g_k * n^4 * (V - V_k) - g_n * n^4 * (V - V_n)
    dK_e = F * g_k * n^4 * (V - V_k) + F * g_n * n^4 * (V - V_n) - gamma_e * (K_e - K_m)
    dK_i = -F * g_k * n^4 * (V - V_n) - F * g_n * n^4 * (V - V_n) 
    # dn = alpha * (1 - n) - beta * n

    dn₀ = -4*α(V)*n₀ + β(V)*n₁
    dn₁ = -(3*α(V) + β(V))*n₁ + 4*α(V)*n₀ + 2*β(V)*n₂
    dn₂ = -(2*α(V) + 2*β(V))*n₂ + 3*α(V)*n₁ + 3*β(V)*n₃
    dn₃ = -(α(V)+3*β(V))*n₃ + 2*α(V)*n₂ + 4*β(V)*n₄
    dn₄ = -4*β(V)*n₄ + α(V)*n₃

    @SVector [dV, dK_e, dK_i, dn₀, dn₁, dn₂, dn₃, dn₄]
end
# Initial conditions
u₀_states=rand(4);

tspan = (0,1000);
h = 1e-3;

sol_states = euler(hodg_hux_det_states, u₀_states, p, tspan, h);

plot(sol_states.t,sol_states.u[1,:],label="V")
plot!(sol_states.t,sol_states.u[2,:],label=L"K_e")
plot!(sol_states.t,sol_states.u[3,:],label=L"K_i")
plot!(sol_states.t,sol_states.u[8,:],label=L"n_4")



# plot(sol_det.t,sol_det.u[1,:],label="V")
# # n
# plot(sol_det.t,sol_det.u[6,:],label=L"n_4") #states
# # m
# plot!(sol_det.t,sol_det.u[10,:],label=L"m_3") #states
# # h
# plot!(sol_det.t,sol_det.u[11,:],label=L"h") #states