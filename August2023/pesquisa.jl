#using .PortfolioOptCTG
using StatsPlots, Plots
using Gurobi

include("C:/Users/frog-/Downloads/Paper/src/Port_Alloc.jl");
include("C:/Users/frog-/Downloads/Paper/src/PortfolioOptCTG.jl");

S = 1200; # Numero de cenários
T = 12; # Numero de meses
tx = 0; # Custo oportunidade anual
tx_m = (1+tx)^(1/12) -1; # Custo oportunidade mensal
#Solv = HiGHS.Optimizer # Solver
Solv = Gurobi.Optimizer # Solver
q = 1 / S .* ones(S) # Probabilidade dos cenarios
pu_Money            = 1e6;

# Leitura do PLD Newave
dir = dirname(dirname(@__FILE__))
PLD = Matrix(DataFrame(CSV.File("PLD.csv";header=false)))'[:,:];
#PLD = Matrix(DataFrame(CSV.File("Monthly_Report/August2023/PLD.csv";header=false)))'[:,:];

#Leitura número de horas no mês
h = Int.(Matrix(DataFrame(CSV.File("horas.csv";header=false)))[:,1]);
#h = Int.(Matrix(DataFrame(CSV.File("Monthly_Report/August2023/horas.csv";header=false)))[:,1]);

usina = ["ARINOS", "DIS"] # Usinas do caso
#usina = ["DIS"]
nI = length(usina) # Número de usinas do caso
gu = zeros(S, T, nI) # Geração das usinas
C = zeros(nI,T)

GF = [50.0, 50.0 # Garantia fisica
]; 

TotGF = sum(GF);

#Leitura geração usina HYDRO (% garantia física)
#gu_dis = Matrix(DataFrame(CSV.File("Monthly_Report/August2023/gu_dis.csv";header=false)))'[1:S,1:T];
gu_dis = Matrix(DataFrame(CSV.File("gu_dis.csv";header=false)))'[1:S,1:T];

gu_arinos = Matrix(DataFrame(CSV.File("gu_arinos.csv";header=false)))'[1:S,1:T];
#gu_arinos = Matrix(DataFrame(CSV.File("Monthly_Report/August2023/gu_arinos.csv";header=false)))'[1:S,1:T];

gu_hydro = Matrix(DataFrame(CSV.File("gu_hydro.csv";header=false)))'[1:S,1:T];
#gu_hydro = Matrix(DataFrame(CSV.File("Monthly_Report/August2023/gu_hydro.csv";header=false)))'[1:S,1:T];

gu = ones(S,T,nI)
gu[:,:,1] = gu_arinos;
gu[:,:,2] = gu_dis;
#gu[:,:,3] = gu_hydro;

ymax = ones(nI); # Porcentagem máxima da geração das usinas contratadas
ymin = ones(nI); # Porcentagem minima da geração das usinas contratada

J = 13 # Número de contratos
Qmax = sum(GF[i] for i = 1:nI);
Q = ones(J,T);
Q[1,:] .= Qmax;
Q[2:13,:] .= LinearAlgebra.I(12).*Qmax;

P = ones(J, T);

a_Y = 0.93;
b_Y = 15.31;

P[1,:] .= a_Y*mean(PLD[1:S,1:T]) + b_Y;

media_m = zeros(T);
media_m = Statistics.mean(PLD, dims = 1);

a_M = a_Y;
b_M = b_Y;

for j in 1:T
    P[j+1,:] .= LinearAlgebra.I(T)[j,:].*(a_M.*media_m[j] + b_M);
end

global gammax = 0.2;
xmax = ones(J)*gammax; # Porcentagem maxima de contratação 
xmin = ones(J).*-gammax;
xmax[1] = 1.0; # Porcentagem minima de contratação
xmin[1] = 0.0; # Porcentagem minima de contratação

α = 0.95 # Percentil CVaR Risco

flag_restricao = false # Boolean restrição receita minima
Rmin = ones(T)*(0.0); # Receita minima

J = 1;

# ================================================================= #
#                          Run Omega Model                          #
# ================================================================= #

#
    # ---> Defining τ/L <---
#

λ = 1.0;

global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
Rdisp_otimo1, Rquant_otimo1, R_otimo1, R_total_otimo1, x_otimo1, y_otimo1 = port_alloc(params);

cvar    = Statistics.mean(sort(R_total_otimo1)[1:Int(round((1 - α) * S)), :]);
avg_adj = mean(R_total_otimo1);
L = cvar/pu_Money;
#L = avg_adj/pu_Money;
h_ = h./pu_Money;

#
    # ---> Solving Proposed Model <---
#

# EXPERIMENTO UPSIDE
model = Model(Solv)

@variable(model, R[1:S])
@variable(model, R_mensal[1:S,1:T])
@variable(model, δ[1:S] >= 0)
@variable(model, x[1:J])
@variable(model, γ >= 0)

@objective(
    model,
    Max,
    sum(q[s]*R[s] for s in 1:S) - L*γ
);

@constraint(model, [s ∈ 1:S], R[s] == sum(R_mensal[s, t] for t in 1:T));

@constraint(model, [s ∈ 1:S], 
    δ[s] >= γ*L - R[s]
);

@constraint(model, [s ∈ 1:S, t ∈ 1:T], 
    R_mensal[s,t] == sum(γ * gu[s, t, i] * GF[i] * PLD[s, t] * h_[t] for i in 1:nI) 
                            + sum((P[j, t] - PLD[s, t]) * Q[j, t] * x[j] * h_[t] for j in 1:J)
);

@constraint(model,[j ∈ 1:J],
    x[j] <= γ*xmax[j]
);

@constraint(model,[j ∈ 1:J],
    x[j] >= γ*xmin[j]
);

@constraint(model,
    sum(q[s] * δ[s] for s in 1:S) == 1
);

@constraint(model, [t ∈ 1:T],
    sum(Q[j, t] * x[j] for j = 1:J) <= sum(GF[i] * (1+gammax) * γ for i = 1:nI)
);

status = optimize!(model);

x_otimo = JuMP.value.(x);
R_otimo = JuMP.value.(R);
γ_otimo = JuMP.value.(γ);
δ_otimo = JuMP.value.(δ);
R_otimo_mensal = JuMP.value.(R_mensal);

x_novo = x_otimo/γ_otimo;
R_otimo_proposto = R_otimo./γ_otimo;
R_otimo_proposto_mensal = R_otimo_mensal./γ_otimo;


#
    # ---> Compute Statistics <---
#

Set_Quant       = [0.01 ; 0.05 ; 0.10 ; 0.25 ; 0.50 ; 0.75 ; 0.90 ; 0.95 ; 0.99];
nQuant          = length(Set_Quant);

global iter_q   = 1;

Quant_Prop      = zeros(1, nQuant);

for q_ ∈ Set_Quant
    Quant_Prop[iter_q] = sort(R_otimo_proposto)[Int(round(q_*S))];
    global iter_q   = iter_q + 1;
end

Avg_Prop        = mean(R_otimo_proposto);
CVaR_Prop       = Statistics.mean(sort(R_otimo_proposto)[1:Int(round((1 - α) * S)), :]);


# ================================================================= #

# ================================================================= #
#                       Run Risk-Reward Model                       #
# ================================================================= #

Set_Λ           = [0.00 ; 0.05 ; 0.10 ; 0.35 ; 0.50 ; 0.75 ; 0.90 ; 0.95 ; 0.99];
n_Λ             = length(Set_Λ);

#Set_Quant       = [0.01 ; 0.05 ; 0.10 ; 0.25 ; 0.50 ; 0.75 ; 0.90 ; 0.95 ; 0.99];
#nQuant          = length(Set_Quant);

global Rdisp_RR        = zeros(n_Λ, S, T, nI);
global Rquant_RR       = zeros(n_Λ, S, T, J);
global RPort_RR        = zeros(n_Λ, S, T);
global RPortTot_RR     = zeros(n_Λ, S);
global xOptimal_RR     = zeros(n_Λ, J);
global yOptimal_RR     = zeros(n_Λ, nI);
global iter_λ          = 1;
global iter_q          = 1;

Quant_RR        = zeros(n_Λ, nQuant);
Avg_RR          = zeros(n_Λ);
CVaR_RR         = zeros(n_Λ);

for λ_ ∈ Set_Λ
    println("\n");
    println("λ = ", λ_);
    println("\n");

    global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
    Rdisp_RR[iter_λ, :, :, :], Rquant_RR[iter_λ, :, :, :], RPort_RR[iter_λ, :, :], RPortTot_RR[iter_λ, :], xOptimal_RR[iter_λ, :], yOptimal_RR[iter_λ, :] = port_alloc(params);

    for q_ ∈ Set_Quant
        Quant_RR[iter_λ, iter_q]    = sort(RPortTot_RR[iter_λ, :])[Int(round(q_*S))]/pu_Money;
        Avg_RR[iter_λ]              = mean(RPortTot_RR[iter_λ, :])/pu_Money;
        CVaR_RR[iter_λ]             = Statistics.mean(sort(RPortTot_RR[iter_λ, :])[1:Int(round((1 - α) * S)), :])/pu_Money;

        global iter_q = iter_q + 1;
    end;

    global iter_q = 1;
    global iter_λ = iter_λ + 1;

end;

# ================================================================= #

# ================================================================= #
#                      Structuring the Results                      #
# ================================================================= #

#
    # ---> Create DataFrames <---
#

xMat    = round.([xOptimal_RR ; x_novo'].*100, digits = 2);
x_df    = DataFrame(xMat, :auto);

println("\n\n");
println(" --> Contracting Levels <-- ");
println("\n");
print(x_df);


qMat    = [Quant_RR ; Quant_Prop];
AvgMat  = [Avg_RR ; Avg_Prop];
CVaRMat = [CVaR_RR ; CVaR_Prop];

qMat    = [AvgMat CVaRMat qMat];
q_df    = DataFrame(qMat, :auto);

println("\n\n");
println(" --> Quantile Levels <-- ");
println("\n");
print(q_df);
println("\n\n\n");


#
    # ---> Plot: Inverse Cumulative Distribution <---
#

s_linewidth         = 2;
pos_Neut            = 1;
pos_Aver            = 4;
λ_plot              = Set_Λ[pos_Aver];

delta = 10;
min_xlim = Int(round(min(minimum(R_otimo_proposto), 
            minimum(RPortTot_RR[pos_Aver, :]./pu_Money), 
            minimum(RPortTot_RR[pos_Neut, :]./pu_Money)
)/10, digits = 0)*10) - delta;

max_xlim = Int(round(max(maximum(R_otimo_proposto), 
            maximum(RPortTot_RR[pos_Aver, :]./pu_Money), 
            maximum(RPortTot_RR[pos_Neut, :]./pu_Money)
)/10, digits = 0)*10) + delta;

p = plot();

R_Prop_Sort         = sort(R_otimo_proposto);

p = plot!((1:S)./S, R_Prop_Sort, ylims = (min_xlim,max_xlim), labels = "Proposed approach",
            legend = :topleft, color=:blue, xlabel="Inverse Cumulative Probability",
            xticks = 0.0:0.1:1.0,
            yticks = min_xlim:10:max_xlim,
            ylabel = "Revenue MMR\$",
            linewidth = s_linewidth
);

R_RR_Sort_Neut      = sort(RPortTot_RR[pos_Neut, :]);

p = plot!((1:S)./S, R_RR_Sort_Neut./pu_Money, 
            color = :black, 
            labels = "Risk-Neutral", 
            linewidth = s_linewidth,
            line=(:dot, 3)
);

R_RR_Sort_Aver      = sort(RPortTot_RR[pos_Aver, :]);

p = plot!((1:S)./S, R_RR_Sort_Aver./pu_Money, 
            color = :gray, 
            labels = string("Risk-Averse (λ = ", λ_plot, ")"),
            linewidth = s_linewidth,
            ls=:dash
);

display(p);
Plots.savefig(p, "InvDistAcum_Full.png");

#
    # ---> Plot: Inverse Cumulative Distribution -- Downside Region <---
#

α_d = 0.20;
Sd = Int(round(S*α_d));

delta = 5;

min_xlim = Int(round(min(minimum(R_Prop_Sort[1:Sd]), 
            minimum(R_RR_Sort_Aver[1:Sd]./pu_Money), 
            minimum(R_RR_Sort_Neut[1:Sd]./pu_Money)
)/10, digits = 0)*10) - delta;

max_xlim = Int(round(max(maximum(R_Prop_Sort[1:Sd]), 
            maximum(R_RR_Sort_Aver[1:Sd]./pu_Money), 
            maximum(R_RR_Sort_Neut[1:Sd]./pu_Money)
)/10, digits = 0)*10) + delta;

p = plot();

R_Prop_Sort         = sort(R_otimo_proposto);

x_axis__            = 0:0.2/239:α_d;

p = plot!(x_axis__, R_Prop_Sort[1:Sd], ylims = (min_xlim,max_xlim), labels = "Proposed approach",
            legend = :bottomright, color=:blue, xlabel="Inverse Cumulative Probability",
            xlims = (0.00,(α_d+0.01)),
            xticks = 0.00:0.02:α_d,
            #yticks = min_xlim:50:max_xlim,
            ylabel = "Revenue MMR\$",
            linewidth = s_linewidth
);

R_RR_Sort_Neut      = sort(RPortTot_RR[pos_Neut, :]);

p = plot!(x_axis__, R_RR_Sort_Neut[1:Sd]./pu_Money, 
            color = :black, 
            labels = "Risk-Neutral", 
            linewidth = s_linewidth,
            line=(:dot, 3)
);

R_RR_Sort_Aver      = sort(RPortTot_RR[pos_Aver, :]);

p = plot!(x_axis__, R_RR_Sort_Aver[1:Sd]./pu_Money, 
            color = :gray, 
            labels = string("Risk-Averse (λ = ", λ_plot, ")"),
            linewidth = s_linewidth,
            ls=:dash
);

display(p);
Plots.savefig(p, "InvDistAcum_DownSide.png");

#
    # ---> Plot: Inverse Cumulative Distribution -- Upside Region <---
#

#α_d = 0.3;
#Sd = Int(round(S*α_d));

p = plot();

delta = 5;

min_xlim = Int(round(min(minimum(R_Prop_Sort[S-Sd+1:S]), 
            minimum(R_RR_Sort_Aver[S-Sd+1:S]./pu_Money), 
            minimum(R_RR_Sort_Neut[S-Sd+1:S]./pu_Money)
)/10, digits = 0)*10) - delta;

max_xlim = Int(round(max(maximum(R_Prop_Sort[S-Sd+1:S]), 
            maximum(R_RR_Sort_Aver[S-Sd+1:S]./pu_Money), 
            maximum(R_RR_Sort_Neut[S-Sd+1:S]./pu_Money)
)/10, digits = 0)*10) + delta;

R_Prop_Sort         = sort(R_otimo_proposto);
x_axis__            = (1 - α_d):0.2/239:1;

#p = plot!((1:Sd)./Sd, R_Prop_Sort[S-Sd+1:S], 
p = plot!(x_axis__, R_Prop_Sort[S-Sd+1:S], 
            ylims = (min_xlim,max_xlim), 
            labels = "Proposed approach",
            legend = :topleft, color=:blue, xlabel="Inverse Cumulative Probability",
            xlims = ((1-α_d),1),
            xticks = (1-α_d):0.02:1,
            yticks = min_xlim:10:max_xlim,
            ylabel = "Revenue MMR\$",
            linewidth = s_linewidth
);

R_RR_Sort_Neut      = sort(RPortTot_RR[pos_Neut, :]);

p = plot!(x_axis__, R_RR_Sort_Neut[S-Sd+1:S]./pu_Money, 
            color = :black, 
            labels = "Risk-Neutral", 
            linewidth = s_linewidth,
            line=(:dot, 3)
);

R_RR_Sort_Aver      = sort(RPortTot_RR[pos_Aver, :]);

p = plot!(x_axis__, R_RR_Sort_Aver[S-Sd+1:S]./pu_Money, 
            color = :gray, 
            labels = string("Risk-Averse (λ = ", λ_plot, ")"),
            linewidth = s_linewidth,
            ls=:dash
);

display(p);
Plots.savefig(p, "InvDistAcum_Upside.png");


#
    # ---> Plot: Contracting Amount Levels <---
#

contrato        = zeros(T,3);

if (J > 1)
    soma_neutral    = xMat[pos_Neut,1] .+ xMat[pos_Neut,2:end];
    soma_averse     = xMat[pos_Aver,1] .+ xMat[pos_Aver,2:end];
    soma_proposto   = xMat[end,1] .+ xMat[end,2:end];
else
    soma_neutral    = xMat[pos_Neut,1];
    soma_averse     = xMat[pos_Aver,1];
    soma_proposto   = xMat[end,1];
end

contrato[:,1]   .= soma_neutral;#*TotGF;
contrato[:,2]   .= soma_averse;#*TotGF;
contrato[:,3]   .= soma_proposto;#*TotGF;

geracao1        = gu_arinos*GF[1];
geracao2        = gu_dis*GF[2];

p_energy = plot();

p_energy = groupedbar!(
    contrato,
    #bar_position = :dodge,
    color = ["black" "gray" "blue"],
    linecolor = :white,
    fillalpha = [0.9 0.9 0.9],
    label = ["Risk-Neutral" string("Risk-Averse (λ = ", λ_plot, ")") "Proposed Approach"],
    #xticks = 1:1:12,
    xlabel = "Bussiness Periods (Months)",
    ylabel = "Energy (avgMW)",
    ylims = (50,300),
    legend = :bottomright,
    #bar_width = 0.80,
    linewidth = 1,
    x_foreground_color_border = :white
);

xaxis_leg       = [i for i ∈ 1:12];
vRES_Generation = zeros(S,T)
for s ∈ 1:S
    for t ∈ 1:T
        vRES_Generation[s,t] = (geracao1[s,t] + geracao2[s,t]);
    end;
end;

lev_fillalpha   = 0.5;
lev_color       = :limegreen;

p_energy.n = 0;

p_energy = violin!(
    vRES_Generation,
    color = lev_color,
    ylims = (0,180),
    legend = :topleft,
    xticks = 1:1:12,
    x_foreground_color_border = :white,
    xlabel = "Bussiness Periods (Months)",
    #ylabel = "Energy Generation (MW)",
    fillalpha=[lev_fillalpha lev_fillalpha lev_fillalpha],
    label = ["Renewable Production" false false false false false false false false false false false false],
    linecolor = nothing
);

display(p_energy);
Plots.savefig(p_energy, "p_energy.png");