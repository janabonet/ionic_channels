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
# const F = 5.6; # membrane capacitance (mM/mV)
F = 5.6/2.7;
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
    alph =  zeros(length(u))
    u[:,1] .= u0
    p[22]=0;
    for i in 1:Int(3600/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
        alph[i] = p[22]
    end
    p[22]=alpha_g;
    # p[22]=4;
    for i in Int(3600/h+1):Int(3780/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
        alph[i] = p[22]
    end
    p[22]=0;
    for i in Int(3780/h+1):Int(10800/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
        alph[i] = p[22]
    end
    p[22]=alpha_g;
    # p[22]=3.4;
    for i in Int(10800/h+1):Int(10980/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
        alph[i] = p[22]
    end
    p[22]=0;
    for i in Int(10980/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
        alph[i] = p[22]
    end
    # p[22]=alpha_g;
    # for i in Int(18000/h+1):Int(18180/h+1)
    #     u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    # end
    # p[22]=0;
    # for i in Int(18180/h+1):n
    #     u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    # end
    return solution_euler(t,u,alph)
end

struct solution_euler
    t::Vector{Float64}
    u::Matrix{Float64}
    alph::Vector{Float64}
end


struct solution_bin
    t::Vector{Float64}
    V::Vector{Float64}
    Ke::Vector{Float64}
    Ki::Vector{Float64}
    N4::Vector{Float64}
end

struct solution_mar
    t::Vector{Float64}
    V::Vector{Float64}
    Ke::Vector{Float64}
    Ki::Vector{Float64}
    N4::Vector{Float64}

    N0::Vector{Float64}
    N1::Vector{Float64}
    N2::Vector{Float64}
    N3::Vector{Float64}
    changes::Vector{Float64}
end
#-------------------------------------------------------------- deterministic, model
function spores_hh_det(u,p,t)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;
    V = u[1]
    K_e = u[2]
    K_i = u[3]
    n = u[4]

    # Nernst potential
    V_k = V_k0*log((K_e/K_i))

    #Dynamical system
    dV = -g_k * n^4 * (V - V_k) - g_n * n^4 * (V - V_n)
    dK_e = F * g_k * n^4 * (V - V_k) + F * g_n * n^4 * (V - V_n) - gamma_e * (K_e - K_m)
    dK_i = -F * g_k * n^4 * (V - V_k) - F * g_n * n^4 * (V - V_n) 
    dn = alpha * (1 - n) - beta * n

    return [dV,dK_e,dK_i,dn]
end

# myrange_det = 1:1000:Int(round(t_tot/h_det));
k_e0=K_m; #0,100,200,400
nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
# k_i0=rand(nd);
k_i0 = 300;
V_0=V_k0*log(k_e0/K_wt) 
#calcul del potencial de nernst per posar-lo com c.i. del potencial de membrana
# u₀_det=[V_0,k_e0,K_wt,rand()];
u₀_det=[V_0wt,k_e0,k_i0,0.0];
tspan = (0,21600);
t_tot=21600;
# t_tot=30;
h_det = 0.5e-3;
myrange_det=1:100:Int(round(t_tot/h_det));
sol_det = euler(spores_hh_det, u₀_det, p, tspan, h_det);

# Plots
f_alph = plot(sol_det.t[myrange_det],sol_det.alph[myrange_det])

f_v=plot(sol_det.t[myrange_det],sol_det.u[1,(myrange_det)], label = "V",xlabel="t (s)", 
ylabel="Membrane potential (mV)",title="k_i0 = "*string(k_i0))

# plot!(xaxis="hores", xticks=0:(h/3600):1000)
fig_conc=plot(sol_det.t[myrange_det],sol_det.u[3,(myrange_det)],label="K_i",xlabel="t (s)", 
ylabel = "Concentration (mM)",title="k_i0 = "*string(k_i0))
plot!(sol_det.t[myrange_det],sol_det.u[2,(myrange_det)],label="K_e",xlabel="t (s)",
ylabel = "Concentration (mM)",title="k_i0 = "*string(k_i0))

# plot(sol_det.t[myrange_det],sol_det.u[4,(myrange_det)],label="n",ylabel="fraction of open subunits")

savefig(f_v,"Vn_ke"*string(k_e0))
savefig(f_v,"Vn_zoom_ke"*string(k_e0))
savefig(fig_conc,"concn_ke"*string(k_e0))
