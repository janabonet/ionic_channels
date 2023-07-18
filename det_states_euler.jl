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
function hodg_hux_det_states(u, p, t)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_tot = p
    # References to variables

    # V = u[1]
    # n = u[2]
    # m = u[3]
    # h = u[4]
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
    # I_na =  g_na * m^3 * h * (V - V_na)
    # I_k  =  g_k * n^4 * (V - V_k)
    # I_l  =  g_l * (V- V_l)
    I_na = g_na * m₃*h * (V - V_na)
    I_k = g_k * n₄ * (V - V_k)
    I_l = g_l * (V - V_l)
   
    # ODE system
    #  dV =  1/C * (I_tot -I_na - I_k - I_l)
    #  dn =  αₙ(V) * (1 - n) - βₙ(V)*n
    #  dm =  αₘ(V) * (1 - m) - βₘ(V)*m
    #  dh =  αₕ(V) * (1 - h) - βₕ(V)*h
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
    return [dV, dn₀, dn₁, dn₂, dn₃, dn₄, dm₀, dm₁, dm₂, dm₃, dh] #11 elements
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

# Initial conditions
n_inf(v) = αₙ(v) / (αₙ(v) + βₙ(v));
m_inf(v) = αₘ(v) / (αₘ(v) + βₘ(v));
h_inf(v) = αₕ(v) / (αₕ(v) + βₕ(v));
v₀ = -60.0;
# u₀ = [v₀, n_inf(v₀), m_inf(v₀), h_inf(v₀)];

n0 = rand(4);
ns0=[n0; n_inf(v₀)];
m0 = rand(3);
ms0 = [m0; m_inf(v₀)]
# h0 = rand();
ns0=ns0/sum(ns0);
ms0=ms0/sum(ms0);
u₀prob = (vcat(v₀,ns0,ms0,h_inf(v₀)));
# u₀prob = SVector{11}(vcat(v₀,ns0,ms0,h_inf(v₀)));
# u₀prob = SVector{11}(vcat(rand(),n0, m0,h0));

# Integration
tspan = (0,200);
h = 1e-3;
sol = euler(hodg_hux_det_states, u₀prob, p, tspan, h);

# Plots
plot(sol.t,sol.u[1,:])
plot(sol.t,sol.u[6,:],label="n")
plot!(sol.t,sol.u[10,:],label="m")
plot!(sol.t,sol.u[11,:],label="h")