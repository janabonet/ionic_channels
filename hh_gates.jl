using Plots, Distributions, LaTeXStrings, BenchmarkTools, DifferentialEquations,StaticArrays
# Parameters
const V_na = 115;
const V_k = -12.0;
const V_l = 10.6;
const g_na = 120;
const g_k = 36.0;
const g_l = 0.3;
const C = 1.0;
I_ext = 0.0;
p = [V_na, V_k, V_l, g_na, g_k, g_l, C, I_ext];

Veq=-0;

αₙ(V) = (0.01 * (10-(V-Veq))) / (exp((10-(V-Veq))/10)-1);
αₘ(V) = (0.1*(25-(V-Veq)))/(exp((25-(V-Veq))/10)-1);
αₕ(V) = 0.07*exp(-(V-Veq)/20);

βₙ(V) = 0.125*exp(-(V-Veq)/80);
βₘ(V) = 4*exp(-(V-Veq)/18);
βₕ(V) = 1/(exp((30-(V-Veq))/10)+1); 

#Time definitions
ti=0;
tf=200;
dt=0.001;
ts=collect(ti:dt:tf);
nt=length(ts);
V=zeros(nt);
m=zeros(nt);
h=zeros(nt);
n=zeros(nt);

# Initial conditions
V[1]= -65; 
m[1]=0.5; 
h[1]=0.06; 
n[1]=0.5; 
I=10.0;

for i in 2:(nt-1)
    V[i+1]=V[i] + dt*(g_na*m[i]^3*h[i]*(V_na-(V[i]))+g_k*n[i]^4*(V_k-V[i])+g_l*(V_l-V[i])+I);
    m[i+1]= m[i] + dt*(αₘ(V[i])*(1-m[i])-βₘ(V[i])*m[i]);
    h[i+1] = h[i] + dt*(αₕ(V[i])*(1-h[i])-βₕ(V[i])*h[i]);
    n[i+1] = n[i] + dt*(αₙ(V[i])*(1-n[i])-βₙ(V[i])*n[i]);
end

fig1=plot(ts,n)
plot!(ts,m)
plot!(ts,h)
