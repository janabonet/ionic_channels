using Plots, Distributions, LaTeXStrings, DifferentialEquations,StaticArrays

# Parameters
const V_na = 115;
const V_k = -12.0;
const V_l = 10.6;
const g_na = 120;
const g_k = 36.0;
const g_l = 0.3;
const C = 1.0;
I_ext = -1;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Gate functions
Veq=-0;
αₙ(V) = (0.01 * (10-(V-Veq))) / (exp((10-(V-Veq))/10)-1);
αₘ(V) = (0.1*(25-(V-Veq)))/(exp((25-(V-Veq))/10)-1);
αₕ(V) = 0.07*exp(-(V-Veq)/20);

βₙ(V) = 0.125*exp(-(V-Veq)/80);
βₘ(V) = 4*exp(-(V-Veq)/18);
βₕ(V) = 1/(exp((30-(V-Veq))/10)+1); 

struct solution_euler
    t::Vector{Float64}
    u::Matrix{Float64}
end

struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N4::Vector{Float64}
    M3::Vector{Float64}
    H::Vector{Float64}
    intensitat::Vector{Float64}
end

struct solution_vars
    t::Vector{Float64}
    V::Vector{Float64}
    N4::Vector{Float64}
    M3::Vector{Float64}
    H::Vector{Float64}

    N0::Vector{Float64}
    N1::Vector{Float64}
    N2::Vector{Float64}
    N3::Vector{Float64}
    changes::Vector{Float64}
    intensitat_vars::Vector{Float64}
end
function euler(f::Function, u0::Vector{Float64}, p::Vector{Float64},
    tspan::Tuple{Int64,Int64}, h::Float64)
    n = round(Int, (tspan[2] - tspan[1]) / h)
    t = collect(tspan[1]:h:tspan[2])
    u = zeros(length(u0), n+1)
    u[:,1] .= u0
    for i in 1:Int(100/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i])
    end
    p[8] += 2.5;
    for i in Int(100/h+1):Int(107/h)
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    p[8] += -2.5;
    for i in Int(107/h+1):n
        u[:,i+1] = u[:,i] + h*f(u[:,i],p,t[i]) #+sqrt(step)*D*rand(d,1) # D es amplitud i d soroll gaussia.
    end
    # p[8] = 0;
    # I_ext=0;
    return solution_euler(t,u)
end

#-------------------------------------------------------------det2
function hodg_hux_det_states(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p
    # References to variables
    V = u[1]

    n₀ = u[2]
    n₁ = u[3]
    n₂ = u[4]
    n₃ = u[5]
    n₄ = u[6]
    m₀ = u[7]
    m₁ = u[8]
    m₂ = u[9]
    m₃ = u[10]
    h = u[11]

    # Channel currents
    I_na = g_na * m₃*h * (V - V_na)
    I_k = g_k * n₄ * (V - V_k)
    I_l = g_l * (V - V_l)
    # I_l=0;

    # ODE system
    dV = 1 / C * (I_ext - I_na - I_k - I_l)

    dn₀ = -4*αₙ(V)*n₀ + βₙ(V)*n₁
    dn₁ = -(3*αₙ(V) + βₙ(V))*n₁ + 4*αₙ(V)*n₀ + 2*βₙ(V)*n₂
    dn₂ = -(2*αₙ(V) + 2*βₙ(V))*n₂ + 3*αₙ(V)*n₁ + 3*βₙ(V)*n₃
    dn₃ = -(αₙ(V)+3*βₙ(V))*n₃ + 2*αₙ(V)*n₂ + 4*βₙ(V)*n₄
    dn₄ = -4*βₙ(V)*n₄ + αₙ(V)*n₃
    
    dm₀ = -3*αₘ(V)*m₀ + βₘ(V)*m₁
    dm₁ = -(2*αₘ(V) + βₘ(V))*m₁ + 3*αₘ(V)*m₀ + 2*βₘ(V)*m₂
    dm₂ = -(αₘ(V) + 2*βₘ(V))*m₂ + 2*αₘ(V)*m₁ + 3*βₘ(V)*m₃
    dm₃ = -3*βₘ(V)*m₃ + αₘ(V)*m₂

    dh = αₕ(V)*(1 - h) - βₕ(V)*h

    @SVector [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀, dm₁, dm₂, dm₃, dh]
end
# -------------------------------------------------- Integration with Euler
# Initial conditions
n_inf(v) = αₙ(v) / (αₙ(v) + βₙ(v));
m_inf(v) = αₘ(v) / (αₘ(v) + βₘ(v));
h_inf(v) = αₕ(v) / (αₕ(v) + βₕ(v));
v₀ = -60.0;
# u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)];

n0 = rand(4);
ns0=[n0; n_inf(v₀)];
m0 = rand(3);
ms0 = [m0; m_inf(v₀)];
# h0 = rand();
ns0=ns0/sum(ns0);
ms0=ms0/sum(ms0);
u₀det = (vcat(v₀,ns0,ms0,h_inf(v₀)));

tspan = (0,300);
h = 1e-3;
sol_det = euler(hodg_hux_det_states, u₀det, p, tspan, h);
# # Plots
# plot(sol_det.t,sol_det.u[1,:],label="V")
# # n
# plot(sol_det.t,sol_det.u[6,:],label=L"n_4") #states
# # m
# plot!(sol_det.t,sol_det.u[10,:],label=L"m_3") #states
# # h
# plot!(sol_det.t,sol_det.u[11,:],label=L"h") #states
# ----------------------------------------------------------------------binomial
I_ext=-4;
function channel_states_bin(N_tot, dt, t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));

    # iniciar vectors
    V = zeros(total_steps)
    intensitat=zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)
    M0 = zeros(total_steps)
    M1 = zeros(total_steps)
    M2 = zeros(total_steps)
    M3 = zeros(total_steps)
    H  = zeros(total_steps)

    # Initial conditions
    V[1] = rand()
    intensitat[1]=0;
    n0 = rand(5)
    m0 = rand(4)
    n0 = round.(n0/sum(n0)*N_tot); 
    m0 = round.(m0/sum(m0)*N_tot); 
    N0[1] = n0[1]
    N1[1] = n0[2]
    N2[1] = n0[3]
    N3[1] = n0[4]
    N4[1] = n0[5]
    M0[1] = m0[1]
    M1[1] = m0[2]
    M2[1] = m0[3]
    M3[1] = m0[4]
    H[1] = round(rand()*N_tot); 

    for i in 2:total_steps
        
        I_ext=-4;
        if i >= 1/dt*100 && i <= 1/dt*107
            I_ext=1.8;
        end

        if i>=1/dt*107
            I_ext=-4;
        end

        intensitat[i]=I_ext;

        # Evolucio canals 
        N0[i] = N0[i-1] + rand(Binomial(N1[i-1],βₙ(V[i-1])*dt)) - rand(Binomial(N0[i-1],4*αₙ(V[i-1])*dt)) 
        N1[i] = N1[i-1] + rand(Binomial(N0[i-1],4*αₙ(V[i-1])*dt)) + rand(Binomial(N2[i-1],2*βₙ(V[i-1])*dt)) - rand(Binomial(N1[i-1],3*αₙ(V[i-1])*dt)) -rand(Binomial(N1[i-1],βₙ(V[i-1])*dt))
        N2[i] = N2[i-1] + rand(Binomial(N1[i-1],3*αₙ(V[i-1])*dt)) + rand(Binomial(N3[i-1],3*βₙ(V[i-1])*dt)) - rand(Binomial(N2[i-1],2*αₙ(V[i-1])*dt)) -rand(Binomial(N2[i-1],2*βₙ(V[i-1])*dt))
        N3[i] = N3[i-1] + rand(Binomial(N2[i-1],2*αₙ(V[i-1])*dt)) + rand(Binomial(N4[i-1],4*βₙ(V[i-1])*dt)) - rand(Binomial(N3[i-1],αₙ(V[i-1])*dt)) - rand(Binomial(N3[i-1],3*βₙ(V[i-1])*dt))
        N4[i] = N4[i-1] + rand(Binomial(N3[i-1],αₙ(V[i-1])*dt)) - rand(Binomial(N4[i-1],4*βₙ(V[i-1])*dt))
       
        M0[i] = M0[i-1] + rand(Binomial(M1[i-1],βₘ(V[i-1])*dt)) - rand(Binomial(M0[i-1],3*αₘ(V[i-1])*dt))
        M1[i] = M1[i-1] + rand(Binomial(M0[i-1],3*αₘ(V[i-1])*dt)) + rand(Binomial(M2[i-1],2*βₘ(V[i-1])*dt)) - rand(Binomial(M1[i-1],2*αₘ(V[i-1])*dt)) - rand(Binomial(M1[i-1],βₘ(V[i-1])*dt))
        M2[i] = M2[i-1] + rand(Binomial(M1[i-1],2*αₘ(V[i-1])*dt)) + rand(Binomial(M3[i-1],3*βₘ(V[i-1])*dt)) - rand(Binomial(M2[i-1],αₘ(V[i-1])*dt)) - rand(Binomial(M2[i-1],2*βₘ(V[i-1])*dt))
        M3[i] = M3[i-1] + rand(Binomial(M2[i-1],αₘ(V[i-1])*dt)) - rand(Binomial(M3[i-1],3*βₘ(V[i-1])*dt))

        H[i] = H[i-1] + rand(Binomial(N_tot .-H[i-1],αₕ(V[i-1])*dt)) - rand(Binomial(H[i-1],βₕ(V[i-1])*dt))

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
        if M0[i] > N_tot
            M0[i]=N_tot
            M1[i]=0;
            M2[i]=0;
            M3[i]=0;
        end
        if M1[i] > N_tot
            M1[i]=N_tot
            M0[i]=0;
            M2[i]=0;
            M3[i]=0;
        end
        if M2[i] > N_tot
            M2[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M3[i]=0;
        end
        if M3[i] > N_tot
            M3[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M2[i]=0;
        end
        if H[i] < 0
            H[i]=0;
            elseif H[i] > N_tot
            H[i]=N_tot
        end
        # Evolucio temporal dels corrents
        I_na = g_na * M3[i-1]/N_tot * H[i-1]/N_tot * (V[i-1] - V_na) ; #println(I_na)
        I_k = g_k * N4[i-1]/N_tot * (V[i-1] - V_k);# println(I_k)
        I_l = g_l * (V[i-1] - V_l); #println(I_l)
        # ODE system
        V[i] = V[i-1] + dt *  1 / C * (I_ext - I_na - I_k - I_l)
    end
    return solution(collect(0:dt:t_tot),V,N4,M3,H,intensitat)
end
#--------------------------------------------------------------------------------------------------------det 2 
I_ext=-4;

function channel_states_markov(N_tot, dt, t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));

    # iniciar vectors
    V = zeros(total_steps)
    intensitat_vars=zeros(total_steps)
    N0 = zeros(total_steps)
    N1 = zeros(total_steps)
    N2 = zeros(total_steps)
    N3 = zeros(total_steps)
    N4 = zeros(total_steps)
    M0 = zeros(total_steps)
    M1 = zeros(total_steps)
    M2 = zeros(total_steps)
    M3 = zeros(total_steps)
    H  = zeros(total_steps)

    changes=zeros(total_steps)

    # Initial conditions
    V[1] = u₀[1];
    intensitat_vars[1]=0;
    # n0 = rand(5)
    n0=u₀[2:6];
    # m0 = rand(4)
    m0=u₀[7:10];
    n0 = round.(n0/sum(n0)*N_tot); 
    m0 = round.(m0/sum(m0)*N_tot); 
    h0=u₀[11];
    N0[1] = n0[1]
    N1[1] = n0[2]
    N2[1] = n0[3]
    N3[1] = n0[4]
    N4[1] = n0[5]
    M0[1] = m0[1]
    M1[1] = m0[2]
    M2[1] = m0[3]
    M3[1] = m0[4]
    H[1] = round(h0*N_tot); 
    
    for i in 2:total_steps

        # t/dt=nº steps = 500/0.5e-5 = 10^8
        # steps/s = 1/dt

        I_ext=-4;
        if i >= 1/dt*100 && i <= 1/dt*107
            I_ext=2.7;
        end

        if i>=1/dt*107
            I_ext=-4;
        end

        intensitat_vars[i]=I_ext;

        #Probabilities definition
        # N0[i] = N0[i-1] +p1c-p0o
        pn1c=βₙ(V[i-1])*dt; 
        pn0o=4*αₙ(V[i-1])*dt;

        # N1[i] = N1[i-1] + p0o + p2c - p1o - p1c
        pn2c = 2*βₙ(V[i-1])*dt;
        pn1o = 3*αₙ(V[i-1])*dt;

        # N2[i] = N2[i-1] + p1o + p3c - p2o - p2c
        pn3c = 3*βₙ(V[i-1])*dt;
        pn2o = 2*αₙ(V[i-1])*dt;

        # N3[i] = N3[i-1] + p2o + p4c - p3o - p3c
        pn4c = 4*βₙ(V[i-1])*dt;
        pn3o = αₙ(V[i-1])*dt;

        # N4[i] = N4[i-1] + p3o - p4c

        # M0[i] = M0[i-1] + pm1c - pm0o
        pm1c = βₘ(V[i-1])*dt;
        pm0o = 3*αₘ(V[i-1])*dt;

        # M1[i] = M1[i-1] + pm0o + pm2c - pm1o - pm1c
        pm2c =2*βₘ(V[i-1])*dt;
        pm1o =2*αₘ(V[i-1])*dt;

        # M2[i] = M2[i-1] + pm1o + pm3c - pm2o - pm2c
        pm3c = 3*βₘ(V[i-1])*dt;
        pm2o = αₘ(V[i-1])*dt;

        # M3[i] = M3[i-1] + p2o - p3c
        
        # H[i] = H[i-1] + pho - phc
        pho = αₕ(V[i-1])*dt;
        phc=βₕ(V[i-1])*dt;
        # phc= (1 - pho)*dt;

        n00=N0[i-1];
        n1=N1[i-1];
        n2=N2[i-1];
        n3=N3[i-1];
        n4=N4[i-1];
        m00=M0[i-1];
        m1=M1[i-1];
        m2=M2[i-1];
        m3=M3[i-1];
        h=H[i-1];

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
        # M channels evolution
        if rand(Uniform(0,1))<pm0o*m00
            m00 = m00 - 1;
            m1 = m1 + 1; 
            ii=ii+1;
        end        
        if rand(Uniform(0,1))<pm1c*m1
            m1 = m1 - 1;
            m00 = m00 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pm1o*m1
            m1 = m1 - 1;
            m2 = m2 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pm2c*m2
            m2 = m2 - 1;
            m1 = m1 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pm2o*m2
            m2 = m2 - 1;
            m3 = m3 + 1;
            ii=ii+1;
        end
        if rand(Uniform(0,1))<pm3c*m3
            m3 = m3 - 1;
            m2 = m2 + 1;
            ii=ii+1;
        end
        # H channels evolution
        if rand(Uniform(0,1))<pho*(N_tot-h)
            h = h + 1;
            ii=ii+1;
        end 
        if rand(Uniform(0,1))< phc*h
            h = h - 1;
            ii=ii+1;
        end
        N0[i]=n00;
        N1[i]=n1;
        N2[i]=n2;
        N3[i]=n3;
        N4[i]=n4;
        M0[i]=m00;
        M1[i]=m1;
        M2[i]=m2;
        M3[i]=m3;
        H[i]=h;
        changes[i-1]=ii;
        # changes[i]=ii;

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
        if M0[i] > N_tot
            M0[i]=N_tot
            M1[i]=0;
            M2[i]=0;
            M3[i]=0;
        end
        if M0[i] < 0
            M0[i]=0
        end
        if M1[i] > N_tot
            M1[i]=N_tot
            M0[i]=0;
            M2[i]=0;
            M3[i]=0;
        end
        if M1[i] < 0
            M1[i]=0
        end
        if M2[i] > N_tot
            M2[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M3[i]=0;
        end
        if M2[i] < 0
            M2[i]=0
        end
        if M3[i] > N_tot
            M3[i]=N_tot
            M0[i]=0;
            M1[i]=0;
            M2[i]=0;
        end
        if M3[i] < 0
            M3[i]=0
        end
        if H[i] < 0
            H[i]=0;
        end
        if H[i] > N_tot
            H[i]=N_tot
        end
        I_na = g_na * M3[i-1]/N_tot * H[i-1]/N_tot * (V[i-1] - V_na) ; #println(I_na)
        I_k = g_k * N4[i-1]/N_tot * (V[i-1] - V_k); 
        I_l = g_l * (V[i-1] - V_l); 
        # I_l=0;

        # ODE system
        V[i] = V[i-1] + dt *  1 / C * (I_ext - I_na - I_k - I_l)
        # println(V[i])
    end
    avg=sum(changes)/total_steps;
    print("mar_avg: "*string(avg))
    return solution_vars(collect(0:dt:t_tot),V,N4,M3,H,N0,N1,N2,N3,changes,intensitat_vars)
end
I_ext=-4;
# -----------------------------------------------------------Simulations
N_tot = 1000;
dt = 0.5e-4;
dt_markov=0.5e-4;
t_tot = 300;
myrange = 1:100:Int(round(t_tot/dt));

# Binomial simulation
sol_bin = channel_states_bin(N_tot, dt, t_tot, p);

# Markov simulation
u₀ = @SVector rand(11);
tspan = (0, 100);
sol_mar = channel_states_markov(N_tot, dt_markov, t_tot, p); #markov solution

# int=plot(sol_mar.t[myrange],sol_mar.intensitat_vars[myrange],label="markov");
# plot!(sol_bin.t[myrange],sol_bin.intensitat[myrange],label="binomial")

#plot vol states
fig1 = plot(sol_mar.t[myrange],sol_mar.V[myrange], label=L"V_{Markov}",
xlabel = L"t (ms)",ylabel = L"V (mV)",dpi=600,size = (700,400))
#plot determinisitc
plot!(sol_det.t, sol_det.u[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)",
linewidth = 1,label=L"V_{det}", ls=:dash,dpi=600,size = (700,400),
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)
# plot binomial
plot!(sol_bin.t[myrange],sol_bin.V[myrange],label=L"V_{bin}",
xlabel = L"t (ms)",ylabel = L"V (mV)",dpi=600,size = (700,400))

fign=plot(sol_mar.t[myrange],sol_mar.N4[myrange],label=L"N_{4,Markov}",dpi=600,size = (700,400))
plot!(sol_det.t,sol_det.u[6,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
linewidth = 1,label=L"n_{det} \cdot N_{tot}",ls=:dash,dpi=600)
plot!(sol_bin.t[myrange],sol_bin.N4[myrange],label=L"N_{4,bin}",dpi=600)

figm=plot(sol_mar.t[myrange],sol_mar.M3[myrange],label=L"M_{3,Markov}",dpi=600,size = (700,400))
plot!(sol_det.t,sol_det.u[10,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
linewidth = 1,label=L"m_{det} \cdot N_{tot}", ls=:dash,dpi=600)
plot!(sol_bin.t[myrange],sol_bin.M3[myrange],label=L"M_{3,bin}",dpi=600)

figh=plot(sol_mar.t[myrange],sol_mar.H[myrange],label=L"H_{Markov}",
xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,size = (700,400),
background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)
plot!(sol_det.t,sol_det.u[11,:]*N_tot,xlabel = "t (ms)", ylabel = "Number of open channels",
linewidth = 1,label=L"h_{det} \cdot N_{tot}", ls=:dash,dpi=600,
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,
legendfontsize=15,left_margin=2Plots.mm, bottom_margin=2Plots.mm)
plot!(sol_bin.t[myrange],sol_bin.H[myrange],label=L"H_{bin}",
xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,
background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)

# fig2 = plot(sol.t[myrange],sol.N4[myrange],label=L"N_{4,Markov}",dpi=600,size = (700,400))
# plot!(sol.t[myrange],sol.M3[myrange],label=L"M_{3,Markov}",dpi=600,size = (700,400))
# plot!(sol.t[myrange],sol.H[myrange],label=L"H_{Markov}",
# xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,size = (700,400),
# background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
# xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)

# #plot states deterministic
# plot!(sol_det.t,sol_det.u[6,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
# linewidth = 1,label=L"n_{det} \cdot N_{tot}",ls=:dash,dpi=600)
# plot!(sol_det.t,sol_det.u[10,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
# linewidth = 1,label=L"m_{det} \cdot N_{tot}", ls=:dash,dpi=600)
# plot!(sol_det.t,sol_det.u[11,:]*N_tot,xlabel = "t (ms)", ylabel = "Number of open channels",
# linewidth = 1,label=L"h_{det} \cdot N_{tot}", ls=:dash,dpi=600,
# xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,
# legendfontsize=15,left_margin=2Plots.mm, bottom_margin=2Plots.mm)

# # plot states binomial
# plot!(sol_bin.t[myrange],sol_bin.N4[myrange],label=L"N_{4,bin}",dpi=600)
# plot!(sol_bin.t[myrange],sol_bin.M3[myrange],label=L"M_{3,bin}",dpi=600)
# plot!(sol_bin.t[myrange],sol_bin.H[myrange],label=L"H_{bin}",
# xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,
# background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
# xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)

# --------------------------------------------Figures
# fig_tot=plot(fig1,fig2,layout=(2,1),dpi=600)
Iext=abs(I_ext);

savefig(fig1,"v_n"*string(N_tot)*"_spike_negin_"*string(Iext))
savefig(fign,"var_n"*string(N_tot)*"_spike_negin_"*string(Iext))
savefig(figm,"var_m"*string(N_tot)*"_spike_negin_"*string(Iext))
savefig(figh,"var_h"*string(N_tot)*"_spike_negin_"*string(Iext))

# savefig(fig1,"v_n"*string(N_tot)*"_spike_negin_2_"*string(Iext))
# savefig(fign,"var_n"*string(N_tot)*"_spike_negin_2_"*string(Iext))
# savefig(figm,"var_m"*string(N_tot)*"_spike_negin_2_"*string(Iext))
# savefig(figh,"var_h"*string(N_tot)*"_spike_negin_2_"*string(Iext))

# savefig(fig2,"nvars_n"*string(N_tot)*"_spike")
a
# for ns in [1000,500,100,50,30,10]
#     sol_n = channel_states_euler(ns, dt, t_tot, p);
    
#     fig1n=plot(sol_n.t[myrange], sol_n.V[myrange], label=L"V_{stoc}",
#     xlabel = L"t (ms)",ylabel = L"V (mV)",dpi=600,size = (700,400))
#     plot!(sol_det.t, sol_det[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)",linewidth = 1,
#     label=L"V_{det}", ls=:dash,dpi=600,size = (700,400),
#     xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)
#     plot!(sol_bin.t[myrange],sol_bin.V[myrange],label=L"V_{bin}",
#     xlabel = L"t (ms)",ylabel = L"V (mV)",dpi=600,size = (700,400))

#     fig2n=plot(sol_n.t[myrange], sol_n.N4[myrange],label=L"N_{4,Markov}",dpi=600,size = (700,400))
#     plot!(sol_n.t[myrange],sol_n.M3[myrange],label=L"M_{3,Markov}",dpi=600,size = (700,400))
#     plot!(sol_n.t[myrange],sol_n.H[myrange],label=L"H_{Markov}",dpi=600,size = (700,400))
    
#     #plot gates deterministic
#     plot!(sol_det.t,sol_det[6,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
#     linewidth = 1,label=L"n_{det} \cdot N_{tot}",ls=:dash,dpi=600)
#     plot!(sol_det.t,sol_det[10,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
#     linewidth = 1,label=L"m_{det} \cdot N_{tot}", ls=:dash,dpi=600)
#     plot!(sol_det.t,sol_det[11,:]*N_tot,xlabel = "t (ms)", ylabel = "Number of open channels",
#     linewidth = 1,label=L"h_{det} \cdot N_{tot}", ls=:dash,dpi=600,
#     xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15,left_margin=2Plots.mm, bottom_margin=2Plots.mm)

#     plot states binomial
#     plot!(sol_bin.t[myrange],sol_bin.N4[myrange],label=L"N_{4,bin}",dpi=600)
#     plot!(sol_bin.t[myrange],sol_bin.M3[myrange],label=L"M_{3,bin}",dpi=600)
#     plot!(sol_bin.t[myrange],sol_bin.H[myrange],label=L"H_{bin}",
#     xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,size = (700,400),
#     background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
#     xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)
# end