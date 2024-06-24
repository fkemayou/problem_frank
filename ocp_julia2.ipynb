using OptimalControl
using CTDirect
using CTBase
using Plots


λ = 300
β =1e-07
μ = 0.05
d = 0.0002
a = 60
r = 0.15
ρ = 400
K = 1000
ν = 0.045
b1 = 7800
b2 = 2
function ocp_T(T, b3)
    @def ocp begin
       t ∈ [ 0, T ], time
       x ∈ R³,  state 
       u ∈ R, control
       S = x₁
       I = x₂
       N = x₃
       ẋ(t) == [ λ-β*(1-u(t))*S(t)*N(t)-μ*S(t), 
                 β*(1-u(t))*S(t)*N(t)-μ*I(t)-d*N(t)*I(t)/(a+I(t)), 
                 (r+ρ*d*I(t)/(a+I(t)))*N(t)*(1-N(t)/(K*I(t)))-ν*N(t) ] # please update! 
       S(0) == 100
       I(0) == 100
       N(0) == 1000
       0 ≤ u(t) ≤ 0.8
       b3*I(T) + ∫(b1*u(t) - b2*S(t) ) → min 
    end
    return ocp
end
# Résolution initiale avec T=15
T_initial = 50
b3_initial = 0
ocp_initial = ocp_T(T_initial, b3_initial)
sol_initial = solve(ocp_initial)
plot(sol_initial)

# 1er Demarrage a chaud
T_new = 100
ocp_new = ocp_T(T_new, b3_initial)
init_new = OCPInit(sol_initial)
sol_new = solve(ocp_new, grid_size = 1000, max_iter = 2000, tol = 10e-20, print_level=5, mu_strategy="adaptive", init=init_new)
plot(sol_new)

# 2e Demarrage a chaud
T_new2 = 200
ocp_new2 = ocp_T(T_new2, b3_initial)
init_new2 = OCPInit(sol_new)
sol_new2 = solve(ocp_new2, grid_size = 1000, max_iter = 2000, tol = 10e-20, print_level=5, mu_strategy="adaptive", init=init_new2)
plot(sol_new2)

# 3e Demarrage a chaud
T_new3 = 300
ocp_new3 = ocp_T(T_new3, b3_initial)
init_new3 = OCPInit(sol_new2)
sol_new3 = solve(ocp_new3, grid_size = 1000, max_iter = 2000, tol = 10e-20, print_level=5, mu_strategy="adaptive", init=init_new3)
plot(sol_new3)

# 4e Demarrage a chaud
T_new4 = 300
b3_new4 = 180
ocp_new4 = ocp_T(T_new4, b3_new4)
init_new4 = OCPInit(sol_new3)
sol_new4 = solve(ocp_new4, grid_size = 1000, max_iter = 2000, tol = 10e-20, print_level=5, mu_strategy="adaptive", init=init_new4)
plot(sol_new4)
