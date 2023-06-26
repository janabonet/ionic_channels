using Plots, Distributions

# HH parameters
const V_na = 55.0;
const V_k = -77.0;
const V_l = -65.0;
const g_na = 40.0;
const g_k = 35.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

# Potassium (n) and sodium (m,h) ion-channel rate functions
function α(V)
    return [
        (0.02 * (V - 25.0)) / (1.0 - exp((-1.0 * (V - 25.0)) / 9.0)),
        (0.182 * (V + 35.0)) / (1.0 - exp((-1.0 * (V + 35.0)) / 9.0)),
        0.25 * exp((-1.0 * (V + 90.0)) / 12.0)
    ]
end
function β(V)
    return [
        (-0.002 * (V - 25.0)) / (1.0 - exp((V - 25.0) / 9.0)),
        (-0.124 * (V + 35.0)) / (1.0 - exp((V + 35.0) / 9.0)),
        (0.25 * exp((V + 62.0) / 6.0)) / exp((V + 90.0) / 12.0)
    ]
end

# Struct to have solution in a concise form
struct solution
    t::Vector{Float64}
    V::Vector{Float64}
    N_o::Matrix{Float64}
    N_c::Matrix{Float64}
end

# Function to simulate 
function hh_gates_simulation(N_tot, dt,t_tot, p)
    V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext = p
    total_steps = Int(t_tot/dt+1);

    # Inicialitzar vectors
    V = zeros(total_steps)
    N_o = zeros(Float64,3,total_steps)
    N_c = zeros(Float64,3,total_steps)

    # Condicions inicials
    V[1] = rand()
    #N_o[:,1] = round.(N_tot*rand(3))

    N_o[:,1] .= 1
    N_c[:,1] .= 0
    #@. N_c[:,1] = N_tot - N_o[:,1]
    # Evolucio temporal
    for i in 2:total_steps
        p1 = β(V[i-1]).*N_o[:,i-1] #probabilitat de transició/unitat de temps (obert-->tancat)
        p2 = α(V[i-1]).*N_c[:,i-1] 

        # a₁ = [rand(Binomial(N_o[:,i-1],p₁_i*dt)) for p₁_i in p₁];
        # a₂ = [rand(Binomial(N_c[:,i-1] ,p₂_i*dt)) for p₂_i in p₂];
        
        a1 = [rand(Binomial(N_o[j,i-1],p1[j]*dt)) for j in 1:length(p1)];
        a2 = [rand(Binomial(N_c[j,i-1] ,p2[j]*dt)) for j in 1:length(p2)];

        # a₁ = rand(Binomial(p₁*dt))
        # a₂ = rand(Binomial(p₂*dt))
        N_o[:,i] = N_o[:,i-1] + (a2-a1); 
        N_c[:,i] = N_c[:,i-1] + (a1-a2); 

        # dno/dt = - alpha no + beta (1-no)
        # no= No/Ntot
        # dNo/dt = - alpha No + beta (Ntot-No)

        # Assegurar que no hi ha estats impossibles
        for j in 1:3
            if N_o[j,i] < 0
                N_o[j,i] = 0 
                N_c[j,i] = N_tot
            elseif N_o[j,i] > N_tot 
                N_o[j,i] = N_tot
                N_c[j,i] = 0
            elseif N_c[j,i] < 0 
                N_o[j,i] = N_tot
                N_c[j,i] = 0
            elseif N_c[j,i] > N_tot 
                N_o[j,i] = 0
                N_c[j,i] = N_tot
            end
        end
        I_na =  g_na * (N_o[2,i-1]/N_tot)^3 * (N_o[3,i-1]/N_tot) * (V[i-1] - V_na)
        I_k  =  g_k * (N_o[1,i-1]/N_tot)^4 * (V[i-1] - V_k)
        I_l  =  g_l * (V[i-1] - V_l)
        V[i] = V[i-1] + dt *(1/C * (I_ext -I_na - I_k - I_l))
    end
    return solution(collect(0:dt:t_tot),V,N_o,N_c)
end

# Parametres simulacio
N_tot = 1; # Nombre total de canals
dt = 0.5e-3;   # time steps
t_tot = 1000; # temps final

# Simulacio
sol = hh_gates_simulation(N_tot, dt,t_tot, p);
myrange = 1:10:Int(round(t_tot/dt));
# Gràfiques
fig1 = plot(sol.t[myrange],sol.V[myrange], label = "V"; xlabel = "t", ylabel = "V")
fig2 = plot(sol.t[myrange],(sol.N_o[1,myrange]/N_tot).^4, label = "n⁴", xlabel = "t")
       plot!(sol.t[myrange],(sol.N_o[2,myrange]/N_tot).^3, label = "m³")
       plot!(sol.t[myrange],sol.N_o[3,myrange]/N_tot, label = "h")

plot(fig1,fig2, layout = (2,1))
