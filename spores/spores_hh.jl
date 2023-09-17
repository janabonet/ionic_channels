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
    p[22]=alpha_g;
    # p[22]=4;
    for i in Int(3600/h+1):Int(3780/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=0;
    for i in Int(3780/h+1):Int(10800/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=alpha_g;
    # p[22]=3.4;
    for i in Int(10800/h+1):Int(10980/h+1)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[22]=0;
    for i in Int(10980/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    # p[22]=alpha_g;
    # for i in Int(18000/h+1):Int(18180/h+1)
    #     u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    # end
    # p[22]=0;
    # for i in Int(18180/h+1):n
    #     u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    # end
    return solution_euler(t,u)
end

struct solution_euler
    t::Vector{Float64}
    u::Matrix{Float64}
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
k_e0=100; #0,100,200,400
nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
k_i0=rand(nd);
V_0=V_k0*log(k_e0/K_wt) 
#calcul del potencial de nernst per posar-lo com c.i. del potencial de membrana
u₀_det=[V_0,k_e0,K_wt,rand()];
# u₀_det=[V_0wt,k_e0,K_wt,rand()];
tspan = (0,21600);
t_tot=21600;
t_tot=30;
h_det = 0.5e-3;
myrange_det=1:100:Int(round(t_tot/h_det));
sol_det = euler(spores_hh_det, u₀_det, p, tspan, h_det);

# Plots
f_v=plot(sol_det.t[myrange_det],sol_det.u[1,(myrange_det)], label = "V",xlabel="t (s)", 
ylabel="Membrane potential (mV)")

# plot!(xaxis="hores", xticks=0:(h/3600):1000)
fig_conc=plot(sol_det.t[myrange_det],sol_det.u[2,(myrange_det)],label="K_e",xlabel="t (s)",
ylabel = "Concentration (mM)")
plot!(sol_det.t[myrange_det],sol_det.u[3,(myrange_det)],label="K_i",xlabel="t (s)", 
ylabel = "Concentration (mM)")
# plot(sol_det.t[myrange_det],sol_det.u[4,(myrange_det)],label="n",ylabel="fraction of open subunits")

savefig(f_v,"Vn_ke"*string(k_e0))
savefig(f_v,"Vn_zoom_ke"*string(k_e0))
savefig(fig_conc,"concn_ke"*string(k_e0))


#---------------------------------------------------------------- deterministic, states
function spores_states(u, p, t)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;
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
    dV = -g_k * n₄^4 * (V - V_k) - g_n * n₄^4 * (V - V_n)
    dK_e = F * g_k * n₄^4 * (V - V_k) + F * g_n * n₄^4 * (V - V_n) - gamma_e * (K_e - K_m)
    dK_i = -F * g_k * n₄^4 * (V - V_n) - F * g_n * n₄^4 * (V - V_n) 
    # dn = alpha * (1 - n) - beta * n

    dn₀ = -4*alpha*n₀ + beta*n₁
    dn₁ = -(3*alpha + beta)*n₁ + 4*alpha*n₀ + 2*beta*n₂
    dn₂ = -(2*alpha + 2*beta)*n₂ + 3*alpha*n₁ + 3*beta*n₃
    dn₃ = -(alpha+3*beta)*n₃ + 2*alpha*n₂ + 4*beta*n₄
    dn₄ = -4*beta*n₄ + alpha*n₃

    @SVector [dV, dK_e, dK_i, dn₀, dn₁, dn₂, dn₃, dn₄]
end
# Initial conditions
n0 = rand(4);
ns0=[n0; rand()];
# h0 = rand();
ns0=ns0/sum(ns0);
k_e0=400;
nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
k_i0=rand(nd);
u₀_states=[V_0wt,k_e0,k_i0,ns0[1],ns0[2],ns0[3],ns0[4],ns0[5]];
tspan = (0,14400);
h = 1e-3;
# t_hores=sol_det.t*h/3600;
# t_hores=t_hores[1:1000];
tot=Int(round(4*3600/h));

sol_states = euler(spores_states, u₀_states, p, tspan, h);

# plot(sol_det.t[1:1000],sol_det.u[1,1:1000], label = "V")
plot(sol_states.t,sol_states.u[1,:],label="V")
plot(sol_states.t,sol_states.u[2,:],label=L"K_e")
plot!(sol_states.t,sol_states.u[3,:],label=L"K_i")
plot(sol_states.t,sol_states.u[8,1:100000],label=L"n_4")
# ---------------------------------------------------------------------------------------- binomial
function spores_states_bin(N_tot, dt, t_tot, p)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));

    # iniciar vectors
    V = zeros(total_steps)
    Ke=zeros(total_steps)
    Ki=zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)

    # Initial conditions
    V[1] = rand()
    k_e0=400;
    Ke[1]=k_e0;
    nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
    k_i0=rand(nd);
    Ki[1]=k_i0;
    n0 = rand(5)
    n0 = round.(n0/sum(n0)*N_tot); 
    N0[1] = n0[1]
    N1[1] = n0[2]
    N2[1] = n0[3]
    N3[1] = n0[4]
    N4[1] = n0[5]

    for i in 2:total_steps

        # Germinant pulses at 1h and 3h
        alpha = 0;
        if (i >=1/dt*3600 && i <= 1/dt*3781) || (i >=1/dt*10801 && i <= 1/dt*10980)
            alpha = alpha_g
        end

        # Evolucio canals 
        N0[i] = N0[i-1] + rand(Binomial(N1[i-1],beta*dt)) - rand(Binomial(N0[i-1],4*alpha*dt)) 
        N1[i] = N1[i-1] + rand(Binomial(N0[i-1],4*alpha*dt)) + rand(Binomial(N2[i-1],2*beta*dt)) - rand(Binomial(N1[i-1],3*alpha*dt)) -rand(Binomial(N1[i-1],beta*dt))
        N2[i] = N2[i-1] + rand(Binomial(N1[i-1],3*alpha*dt)) + rand(Binomial(N3[i-1],3*beta*dt)) - rand(Binomial(N2[i-1],2*alpha*dt)) -rand(Binomial(N2[i-1],2*beta*dt))
        N3[i] = N3[i-1] + rand(Binomial(N2[i-1],2*alpha*dt)) + rand(Binomial(N4[i-1],4*beta*dt)) - rand(Binomial(N3[i-1],alpha*dt)) - rand(Binomial(N3[i-1],3*beta*dt))
        N4[i] = N4[i-1] + rand(Binomial(N3[i-1],alpha*dt)) - rand(Binomial(N4[i-1],4*beta*dt))
       
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

        V_k = V_k0*log(abs(K_e/K_i))
        # ODE system
        V[i] = V[i-1] + dt * (-g_k * N4[i-1]^4 * (V[i-1] - V_k) - g_n * N4[i-1]^4 * (V[i-1] - V_n))
        Ke[i] = Ke[i-1] + dt * (F * g_k * N4[i-1]^4 * (V[i-1] - V_k) + F * g_n * N4[i-1]^4 * (V[i-1] - V_n) - gamma_e * (Ke[i-1] - K_m))
        Ki[i] = Ki[i-1] + dt * (-F * g_k * N4[i-1]^4 * (V[i-1] - V_n) - F * g_n * N4[i-1]^4 * (V[i-1] - V_n) )
    end
    return solution_bin(collect(0:dt:t_tot),V,Ke,Ki,N4)
end
#--------------------------------------------------------------------------------------------------------det 2 
alpha=0;
function channel_states_markov(N_tot, dt, t_tot, p)
    g_k, g_kq, g_n, g_nq, V_k0, V_n, alpha_g, beta, V_0wt, V_0ktrc,
    V_0yugO, gamma_e, F, K_m, K_s, K_wt, K_ktrc, K_yugO, sigma_wt, sigma_ktrc, sigma_yugO,alpha = p;

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));

    # iniciar vectors
    V = zeros(total_steps)
    Ke=zeros(total_steps)
    Ki=zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)
    changes=zeros(total_steps)

    # Initial conditions
    V[1] = rand()
    k_e0=400;
    Ke[1]=k_e0;
    nd=Normal(K_wt,sigma_wt); # Distribució normal per la concentració intracelular
    k_i0=rand(nd);
    Ki[1]=k_i0;
    n0 = rand(5)
    n0 = round.(n0/sum(n0)*N_tot); 
    N0[1] = n0[1]
    N1[1] = n0[2]
    N2[1] = n0[3]
    N3[1] = n0[4]
    N4[1] = n0[5]
    
    for i in 2:total_steps

        # Germinant pulses at 1h and 3h
        alpha = 0;
        if (i >=1/dt*3600 && i <= 1/dt*3781) || (i >=1/dt*10801 && i <= 1/dt*10980)
            alpha = alpha_g
        end

        #Probabilities definition
        # N0[i] = N0[i-1] +p1c-p0o
        pn1c=beta*dt; 
        pn0o=4*alpha*dt;

        # N1[i] = N1[i-1] + p0o + p2c - p1o - p1c
        pn2c = 2*beta*dt;
        pn1o = 3*alpha*dt;

        # N2[i] = N2[i-1] + p1o + p3c - p2o - p2c
        pn3c = 3*beta*dt;
        pn2o = 2*alpha*dt;

        # N3[i] = N3[i-1] + p2o + p4c - p3o - p3c
        pn4c = 4*beta*dt;
        pn3o = alpha*dt;

        # N4[i] = N4[i-1] + p3o - p4c

        n00=N0[i-1];
        n1=N1[i-1];
        n2=N2[i-1];
        n3=N3[i-1];
        n4=N4[i-1];
        ii=changes[i-1];
        
        # N channels evolution
        if rand(Uniform(0,1))<pn0o*n00
            # n0 = N0[i-1] - 1;
            n0 = n00 - 1;
            n1 = n1 + 1; 
            ii=ii+1;
        end        
        if rand(Uniform(0,1))<pn1c*n1
            n1 = n1- 1;
            n0 = n00 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn1o*n1
            n1 = n1 - 1;
            n2 = n2 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn2c*n2
            n2 = n2 - 1;
            n1 = n1 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn2o*n2
            n2 = n2 - 1;
            n3 = n3 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn3c*n3
            n3 = n3 - 1;
            n2 = n2 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn3o*n3
            n3 = n3 - 1;
            n4 = n4 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pn4c*n4
            n4 = n4 - 1;
            n3 = n3 + 1;
            ii=ii+1;
        end
        # Values update
        N0[i]=n00;
        N1[i]=n1;
        N2[i]=n2;
        N3[i]=n3;
        N4[i]=n4;
        changes[i-1]=ii;

        # Condition to avoid physically impossible states
        if N0[i] > N_tot
            N0[i]=N_tot
            N1[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N0[i] < 0
            N0[i]=0
        end
        if N1[i] > N_tot
            N1[i]=N_tot
            N0[i]=0;
            N2[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N1[i] < 0
            N1[i]=0
        end
        if N2[i]>N_tot
            N2[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N3[i]=0;
            N4[i]=0;
        end
        if N2[i] < 0
            N2[i]=0
        end
        if N3[i] > N_tot
            N3[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end
        if N3[i] < 0
            N3[i]=0
        end
        if N4[i] > N_tot
            N4[i]=N_tot
            N0[i]=0;
            N1[i]=0;
            N2[i]=0;
        end
        if N4[i] < 0
            N4[i]=0
        end

        V_k = V_k0*log(abs(K_e/K_i))
        # ODE system
        V[i] = V[i-1] + dt * (-g_k * N4[i-1]^4 * (V[i-1] - V_k) - g_n * N4[i-1]^4 * (V[i-1] - V_n))
        Ke[i] = Ke[i-1] + dt * (F * g_k * N4[i-1]^4 * (V[i-1] - V_k) + F * g_n * N4[i-1]^4 * (V[i-1] - V_n) - gamma_e * (Ke[i-1] - K_m))
        Ki[i] = Ki[i-1] + dt * (-F * g_k * N4[i-1]^4 * (V[i-1] - V_n) - F * g_n * N4[i-1]^4 * (V[i-1] - V_n) )
    end
    avg=sum(changes)/total_steps;
    print("mar_avg: "*string(avg))
    return solution_mar(collect(0:dt:t_tot),V,Ke,Ki,N4,N0,N1,N2,N3,changes)
end
alpha=0;
# -----------------------------------------------------------Simulations
N_tot = 500;
dt = 0.5e-4;
dt_markov=0.5e-4;
t_tot = 14400;
myrange = 1:100:Int(round(t_tot/dt));

# Binomial simulation
sol_bin = spores_states_bin(N_tot, dt, t_tot, p);

# Markov simulation
u₀ = @SVector rand(11);
tspan = (0, 14400);
sol_mar = spores_states_markov(N_tot, dt_markov, t_tot, p); #markov solution