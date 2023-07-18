using Plots, Distributions, LaTeXStrings, DifferentialEquations,StaticArrays

# Integration functions and structs
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

# DIff eqs system
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
     dn =  αₙ(V) * (1 - n) - βₙ(V)*n
     dm =  αₘ(V) * (1 - m) - βₘ(V)*m
     dh =  αₕ(V) * (1 - h) - βₕ(V)*h
    return [dV,dn,dm,dh]
end 

# Parameters
const V_na = 115;
const V_k = -12.0;
const V_l = 10.6;
const g_na = 120;
const g_k = 36.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 9.8;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

Veq=-0;

αₙ(V) = (0.01 * (10-(V-Veq))) / (exp((10-(V-Veq))/10)-1);
αₘ(V) = (0.1*(25-(V-Veq)))/(exp((25-(V-Veq))/10)-1);
αₕ(V) = 0.07*exp(-(V-Veq)/20);

βₙ(V) = 0.125*exp(-(V-Veq)/80);
βₘ(V) = 4*exp(-(V-Veq)/18);
βₕ(V) = 1/(exp((30-(V-Veq))/10)+1); 

# ----------------------------------Integration with loop euler method

#Time definitions
# ti=0;
# tf=200;
# dt=0.001;
# ts=collect(ti:dt:tf);
# nt=length(ts);
# V=zeros(nt);
# m=zeros(nt);
# h=zeros(nt);
# n=zeros(nt);

# V[1]= -65; 
# m[1]=0.5; 
# h[1]=0.06; 
# n[1]=0.5; 
# I=10.0;
# for i in 2:(nt-1)
#     V[i+1]=V[i] + dt*(g_na*m[i]^3*h[i]*(V_na-(V[i]))+g_k*n[i]^4*(V_k-V[i])+g_l*(V_l-V[i])+I);
#     m[i+1]= m[i] + dt*(αₘ(V[i])*(1-m[i])-βₘ(V[i])*m[i]);
#     h[i+1] = h[i] + dt*(αₕ(V[i])*(1-h[i])-βₕ(V[i])*h[i]);
#     n[i+1] = n[i] + dt*(αₙ(V[i])*(1-n[i])-βₙ(V[i])*n[i]);
# end

# fig1=plot(ts,n)
# plot!(ts,m)
# plot!(ts,h)

# ----------------------------------------------------------------------
# Initial conditions
n_inf(v) = αₙ(v) / (αₙ(v) + βₙ(v));
m_inf(v) = αₘ(v) / (αₘ(v) + βₘ(v));
h_inf(v) = αₕ(v) / (αₕ(v) + βₕ(v));
v₀ = -60.0;
u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)];

tspan = (0,200);
h = 1e-3;
sol = euler(hodg_hux_det, u₀, p, tspan, h);

# Plots
plot(sol.t,sol.u[1,:])
plot(sol.t,sol.u[2,:],label="n")
plot!(sol.t,sol.u[3,:],label="m")
plot!(sol.t,sol.u[4,:],label="h")