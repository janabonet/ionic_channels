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
# F = 5.6/2.5;
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

pols=3; #durada del pols (s'ha de canviar també dins la funció)

struct solution_bin
    t::Vector{Float64}
    V::Vector{Float64}
    Ke::Vector{Float64}
    Ki::Vector{Float64}
    alphas::Vector{Float64}
    N0::Vector{Float64}
    N1::Vector{Float64}
    N2::Vector{Float64}
    N3::Vector{Float64}
    N4::Vector{Float64}
end

# ---------------------------------------------------------------------------------------- binomial
function spores_states_bin(N_tot, dt, t_tot, p)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));
    # total_steps = 2000;
    # iniciar vectors
    V = zeros(total_steps)
    Ke=zeros(total_steps)
    Ki=zeros(total_steps)
    alphas = zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)
    # Initial conditions
    V[1] = V_0wt;
    k_e0=K_m;
    K_e = K_m;
    Ke[1]=k_e0;
    nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
    # k_i0=rand(nd);
    k_i0 = 300;
    Ki[1]=k_i0; 
    N0[1] = 1*N_tot;
    N1[1] = 0;
    N2[1] = 0;
    N3[1] = 0;
    N4[1] = 0;
    pols = 3; # durada del pols

    for i in 2:total_steps
        # Germinant pulses at 1h and 3h
        alpha = 0;
        if (i >=1/dt*60 && i <= 1/dt*(60+pols)) || (i >=1/dt*180 && i <= 1/dt*(180+pols))|| (i >=1/dt*300 && i <= 1/dt*(300+pols))
            alpha = alpha_g
            alphas[i] = alpha_g;
        end
        # Evolucio canals 
        p0o = rand(Binomial(N0[i-1],4*alpha*dt))
        p1c = rand(Binomial(N1[i-1],beta*dt))
        p1o = rand(Binomial(N1[i-1],3*alpha*dt))
        p2c = rand(Binomial(N2[i-1],2*beta*dt))
        p2o =rand(Binomial(N2[i-1],2*alpha*dt))
        p3c =rand(Binomial(N3[i-1],3*beta*dt))
        p3o =rand(Binomial(N3[i-1],alpha*dt))
        p4c =rand(Binomial(N4[i-1],4*beta*dt))
        N0[i] = N0[i-1] + p1c - p0o; 
        N1[i] = N1[i-1] + p0o + p2c - p1o -p1c
        N2[i] = N2[i-1] + p1o + p3c - p2o - p2c
        N3[i] = N3[i-1] + p2o + p4c - p3o - p3c
        N4[i] = N4[i-1] + p3o - p4c
       
        # Evitem que hi hagi estats sense sentit físic
        if N0[i] > N_tot
            N0[i]=N_tot
            N1[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N1[i] > N_tot
            N1[i]=N_tot
            N0[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N2[i]>N_tot
            N2[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N3[i] > N_tot
            N3[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end
        if N4[i] > N_tot
            N4[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end

        K_i = Ki[i-1];
        V_k = V_k0*log(K_e/K_i)
        # ODE system
        V[i] = V[i-1] + dt * (-g_k * N4[i-1] * (V[i-1] - V_k) - g_n * N4[i-1] * (V[i-1] - V_n))
        Ke[i] = Ke[i-1] + dt * (F * g_k * N4[i-1] * (V[i-1] - V_k) + F * g_n * N4[i-1] * (V[i-1] - V_n) - gamma_e * (Ke[i-1] - K_m))
        Ki[i] = Ki[i-1] + dt * (-F * g_k * N4[i-1] * (V[i-1] - V_k) - F * g_n * N4[i-1] * (V[i-1] - V_n) )
    end
    return solution_bin(collect(0:dt:t_tot),V,Ke,Ki,alphas,N0,N1,N2,N3,N4)
end

# -----------------------------------------------------------Simulations
N_tot = 8;  # número de canals simulats
dt = 0.5e-4;
t_tot = 360;
myrange_bin = 1:1000:Int(round(t_tot/dt));

# Binomial simulation
sol_bin = spores_states_bin(N_tot, dt, t_tot, p);
# myrange_bin = 1:1000:length(sol_bin.t);

# Plots
k_i0 = 300;

fig_alpha = plot(sol_bin.t[myrange_bin],sol_bin.alphas[myrange_bin])

f_v=plot(sol_bin.t[myrange_bin],sol_bin.V[myrange_bin], label = L"V",xlabel="t (s)", 
ylabel=L"Membrane\; potential (mV)",title="k_i0 = "*string(k_i0))


fig_conc=plot(sol_bin.t[myrange_bin],sol_bin.Ki[myrange_bin],label=L"K_i",xlabel=L"t\;(s)", 
ylabel = L"Concentration\;(mM)",title=L"k_i0 = "*string(k_i0)*"_polsos = "*string(pols)*" s")
plot!(sol_bin.t[myrange_bin],sol_bin.Ke[myrange_bin],label=L"K_e",xlabel=L"t\;(s)",
ylabel = L"Concentration \;(mM)",title="k_i0 = "*string(k_i0))

# fig ns
fig_ns=plot(sol_bin.t[myrange_bin],sol_bin.N0[myrange_bin],label=L"N0",
xlabel="t (s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N1[myrange_bin],label=L"N1",
xlabel="t (s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N2[myrange_bin],label=L"N2",
xlabel="t (s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N3[myrange_bin],label=L"N3",
xlabel="t (s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N4[myrange_bin],label=L"N4",
xlabel=L"t \;(s)",title="k_i0 = "*string(k_i0)*", polsos = "*string(pols)*" s",
ylabel=L"Number\;of \;open \;channels",xlimits=(50,80))

fig_nszoom= plot(sol_bin.t[myrange_bin],sol_bin.N1[myrange_bin],label=L"N1",
xlabel=L"t\; (s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N2[myrange_bin],label=L"N2",
xlabel=L"t \;(s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N3[myrange_bin],label=L"N3",
xlabel=L"t \;(s)",title="Ns")
plot!(sol_bin.t[myrange_bin],sol_bin.N4[myrange_bin],label=L"N4",
xlabel=L"t \;(s)",title="k_i0 = "*string(k_i0)*", polsos = "*string(pols)*" s",
xlimits=(50,80),ylabel=L"Number\;of \;open \;channels")


savefig(f_v,"v"*string(N_tot)*"_bin_ki"*string(k_i0)*"_polsos"*string(pols))

savefig(fig_conc,"conc"*string(N_tot)*"_bin_ki"*string(k_i0)*"_polsos"*string(pols))
# savefig(fig_conczoom,"conczoom"*string(N_tot)*"_bin_ki"*string(k_i0)*"_polsos"*string(pols))


savefig(fig_ns,"ns"*string(N_tot)*"_bin_ki"*string(k_i0)*"_polsos"*string(pols))
savefig(fig_nszoom,"nszoom"*string(N_tot)*"_bin_ki"*string(k_i0)*"_polsos"*string(pols))
