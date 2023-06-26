using Plots, Distributions, LaTeXStrings, BenchmarkTools, DifferentialEquations,StaticArrays
# Parameters
const V_na = 55.0;
const V_k = -77.0;
const V_l = -65.0;
const g_na = 40.0;
const g_k = 35.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];
# Gate functions
αₙ(V) = (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0));
αₘ(V) = (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0));
αₕ(V) = 0.25 * exp((-1.0 * (V + 90.0)) / 12.0);

βₙ(V) = (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0));
βₘ(V) = (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0));
βₕ(V) = (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0);


struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N4::Vector{Float64}
    M3::Vector{Float64}
    H::Vector{Float64}
end

function channel_states_euler(N_tot, dt, t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p

    # integration Parameters
    total_steps = Int(round(t_tot/dt+1));

    # iniciar vectors
    V = zeros(total_steps)
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
        #Probabilities definition
        #N0[i] = N0[i-1] + rand(Binomial(N1[i-1],βₙ(V[i-1])*dt)) - rand(Binomial(N0[i-1],4*αₙ(V[i-1])*dt)) 
        # N0[i] = N0[i-1] +p1c-p0o
        pn1c=βₙ(V[i-1])*dt; 
        pn0o=4*αₙ(V[i-1])*dt;
        # N1[i] = N1[i-1] + rand(Binomial(N0[i-1],4*αₙ(V[i-1])*dt)) + rand(Binomial(N2[i-1],2*βₙ(V[i-1])*dt)) - rand(Binomial(N1[i-1],3*αₙ(V[i-1])*dt)) -rand(Binomial(N1[i-1],βₙ(V[i-1])*dt))
        # N1[i] = N1[i-1] + p0o + p2c - p1o - p1c
        pn2c = 2*βₙ(V[i-1])*dt;
        pn1o = 3*αₙ(V[i-1])*dt;
        # N2[i] = N2[i-1] + rand(Binomial(N1[i-1],3*αₙ(V[i-1])*dt)) + rand(Binomial(N3[i-1],3*βₙ(V[i-1])*dt)) - rand(Binomial(N2[i-1],2*αₙ(V[i-1])*dt)) -rand(Binomial(N2[i-1],2*βₙ(V[i-1])*dt))
        # N2[i] = N2[i-1] + p1o + p3c - p2o - p2c
        pn3c = 3*βₙ(V[i-1])*dt;
        pn2o = 2*αₙ(V[i-1])*dt;
        # N3[i] = N3[i-1] + rand(Binomial(N2[i-1],2*αₙ(V[i-1])*dt)) + rand(Binomial(N4[i-1],4*βₙ(V[i-1])*dt)) - rand(Binomial(N3[i-1],αₙ(V[i-1])*dt)) - rand(Binomial(N3[i-1],3*βₙ(V[i-1])*dt))
        # N3[i] = N3[i-1] + p2o + p4c - p3o - p3c
        pn4c = 4*βₙ(V[i-1])*dt;
        pn3o = αₙ(V[i-1])*dt;
        # N4[i] = N4[i-1] + rand(Binomial(N3[i-1],αₙ(V[i-1])*dt)) - rand(Binomial(N4[i-1],4*βₙ(V[i-1])*dt))
        # N4[i] = N4[i-1] + p3o - p4c

        # M0[i] = M0[i-1] + rand(Binomial(M1[i-1],βₘ(V[i-1])*dt)) - rand(Binomial(M0[i-1],3*αₘ(V[i-1])*dt))
        # M0[i] = M0[i-1] + pm1c - pm0o
        pm1c = βₘ(V[i-1])*dt;
        pm0o = 3*αₘ(V[i-1])*dt;
        # M1[i] = M1[i-1] + rand(Binomial(M0[i-1],3*αₘ(V[i-1])*dt)) + rand(Binomial(M2[i-1],2*βₘ(V[i-1])*dt)) - rand(Binomial(M1[i-1],2*αₘ(V[i-1])*dt)) + rand(Binomial(M1[i-1],βₘ(V[i-1])*dt))
        # M1[i] = M1[i-1] + pm0o + pm2c - pm1o - pm1c
        pm2c =2*βₘ(V[i-1])*dt;
        pm1o =2*αₘ(V[i-1])*dt;
        # M2[i] = M2[i-1] + rand(Binomial(M1[i-1],2*αₘ(V[i-1])*dt)) + rand(Binomial(M3[i-1],3*βₘ(V[i-1])*dt)) - rand(Binomial(M2[i-1],αₘ(V[i-1])*dt)) - rand(Binomial(M2[i-1],2*βₘ(V[i-1])*dt))
        # M2[i] = M2[i-1] + pm1o + pm3c - pm2o - pm2c
        pm3c = 3*βₘ(V[i-1])*dt;
        pm2o = αₘ(V[i-1])*dt;
        # M3[i] = M3[i-1] + rand(Binomial(M2[i-1],αₘ(V[i-1])*dt)) - rand(Binomial(M3[i-1],3*βₘ(V[i-1])*dt))
        # M3[i] = M3[i-1] + p2o - p3c
        
        # H[i] = H[i-1] + rand(Binomial(N_tot .-H[i-1],αₕ(V[i-1])*dt)) - rand(Binomial(H[i-1],βₕ(V[i-1])*dt))
        # H[i] = H[i-1] + pho - phc
        pho = αₕ(V[i-1])*dt;
        # phc= 1 - pho;
        phc = βₕ(V[i-1])*dt;

 
        # Evolucio canals 
        if rand()<pn0o
            N0[i] = N0[i-1] - 1;
            N1[i] = N1[i-1] + 1; #?? aquests no hi són al seu codi
        end        
        if rand()<pn1c 
            N1[i] = N1[i-1]- 1;
            N0[i] = N0[i-1] + 1;
        end
        if rand()<pn1o
            N1[i] = N1[i-1] - 1;
            N2[i] = N2[i-1] + 1;
        end
        if rand()<pn2c
            N2[i] = N2[i-1] - 1;
            N1[i] = N1[i-1] + 1;
        end
        if rand()<pn2o
            N2[i] = N2[i-1] - 1;
            N3[i] = N3[i-1] + 1;
        end
        if rand()<pn3c
            N3[i] = N3[i-1] - 1;
            N2[i] = N2[i-1] + 1;
        end
        if rand()<pn3o
            N3[i] = N3[i-1] - 1;
            N4[i] = N4[i-1] + 1;
        end

        if rand()<pn4c
            N4[i] = N4[i-1] - 1;
            N3[i] = N3[i-1] + 1;
        end

        # # Evitem que hi hagi estats sense sentit físic
        # if N0[i] > N_tot
        #     N0[i]=N_tot
        #     N1[i]=0;
        #     N2[i]=0;
        #     N3[i]=0;
        #     N4[i]=0;
        #     # elseif N0[i] < 0
        #     #     N0[i]=0;
        # end
        # if N1[i] > N_tot
        #     N1[i]=N_tot
        #     N0[i]=0;
        #     N2[i]=0;
        #     N3[i]=0;
        #     N4[i]=0;
        # end
        # if N2[i]>N_tot
        #     N2[i]=N_tot
        #     N0[i]=0;
        #     N1[i]=0;
        #     N3[i]=0;
        #     N4[i]=0;
        # end
        # if N3[i] > N_tot
        #     N3[i]=N_tot
        #     N0[i]=0;
        #     N1[i]=0;
        #     N2[i]=0;
        # end
        # if N4[i] > N_tot
        #     N4[i]=N_tot
        #     N0[i]=0;
        #     N1[i]=0;
        #     N2[i]=0;
        # end
        # if M0[i] > N_tot
        #     M0[i]=N_tot
        #     M1[i]=0;
        #     M2[i]=0;
        #     M3[i]=0;
        # end
        # if M1[i] > N_tot
        #     M1[i]=N_tot
        #     M0[i]=0;
        #     M2[i]=0;
        #     M3[i]=0;
        # end
        # if M2[i] > N_tot
        #     M2[i]=N_tot
        #     M0[i]=0;
        #     M1[i]=0;
        #     M3[i]=0;
        # end
        # if M3[i] > N_tot
        #     M3[i]=N_tot
        #     M0[i]=0;
        #     M1[i]=0;
        #     M2[i]=0;
        # end
        # if H[i] < 0
        #     H[i]=0;
        #     elseif H[i] > N_tot
        #     H[i]=N_tot
        # end

        I_na = g_na * M3[i-1]/N_tot * H[i-1]/N_tot * (V[i-1] - V_na) ; #println(I_na)
        I_k = g_k * N4[i-1]/N_tot * (V[i-1] - V_k);# println(I_k)
        I_l = g_l * (V[i-1] - V_l); #println(I_l)

        # ODE system
        V[i] = V[i-1] + dt *  1 / C * (I_ext - I_na - I_k - I_l)
        # println(V[i])
    end

    return solution(collect(0:dt:t_tot),V,N4,M3,H)
    
end

#-------------------------------------------------------------det2
function hodg_hux_gates(u, p, t)
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

p[8] = 0.0;

u₀ = @SVector rand(11);
n0 = rand(5);
m0 = rand(4);
h0 = rand();
n0=n0/sum(n0);
m0=m0/sum(m0);
u₀ = SVector{11}(vcat(rand(),n0, m0,h0));
tspan = (0, 500);

# Integration (states)
step_current= PresetTimeCallback(100,integrator -> integrator.p[8] += 0);
prob_det = ODEProblem(hodg_hux_gates, u₀, tspan, p, dtmax = 0.001);
sol_det = solve(prob_det, saveat = 0.1, callback = step_current);
p[8] = 0.0

#--------------------------------------------------------------------------------------------------------det 2 


N_tot = 100;
dt = 0.5e-3;
t_tot = 500;

sol = channel_states_euler(N_tot, dt, t_tot, p);
myrange = 1:100:Int(round(t_tot/dt));

#plot vol states
fig1 = plot(sol.t[myrange],sol.V[myrange], label=L"V_{stoc}",
xlabel = L"t (ms)",ylabel = L"V (mV)",dpi=600,size = (700,400))
#plot determinisitc
plot!(sol_det.t, sol_det[1,:], xlabel = L"t (ms)", ylabel = L"V (mV)",
linewidth = 1,label=L"V_{det}", ls=:dash,dpi=600,size = (700,400),
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)


fig2 = plot(sol.t[myrange],sol.N4[myrange],label=L"N_{4,Markov}",dpi=600,size = (700,400))
plot!(sol.t[myrange],sol.M3[myrange],label=L"M_{3,Markov}",dpi=600,size = (700,400))
plot!(sol.t[myrange],sol.H[myrange],label=L"H_{Markov}",
xlabel =L"t (ms)", ylabel ="Number of open channels",dpi=600,size = (700,400),
background_color_legend = :white, foreground_color_legend = nothing,legend=:topright,
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15)

#plot gates deterministic
plot!(sol_det.t,sol_det[6,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
linewidth = 1,label=L"n_{det} \cdot N_{tot}",ls=:dash,dpi=600)
plot!(sol_det.t,sol_det[10,:]*N_tot,xlabel = L"t (ms)", ylabel = L"Number\:of\:open\:channels",
linewidth = 1,label=L"m_{det} \cdot N_{tot}", ls=:dash,dpi=600)
plot!(sol_det.t,sol_det[11,:]*N_tot,xlabel = "t (ms)", ylabel = "Number of open channels",
linewidth = 1,label=L"h_{det} \cdot N_{tot}", ls=:dash,dpi=600,
xtickfontsize=12,ytickfontsize=12,xguidefontsize=16,yguidefontsize=16,legendfontsize=15,left_margin=2Plots.mm, bottom_margin=2Plots.mm)


# fig_tot=plot(fig1,fig2,layout=(2,1),dpi=600)
# savefig(fig1,"v_n"*string(N_tot))
# savefig(fig2,"vars_n"*string(N_tot))