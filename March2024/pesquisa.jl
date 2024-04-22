using StatsPlots, Plots
using HiGHS, LinearAlgebra, XLSX

include("src/Port_Alloc.jl");
include("src/PortfolioOptCTG.jl");

# ------------------------ PARÂMETROS DO CASO ------------------------

S         = 1200;               # Numero de cenários
T         = 12;                 # Numero de meses
tx        = 0;                  # Custo oportunidade anual
tx_m      = (1+tx)^(1/12) -1;   # Custo oportunidade mensal
Solv      = HiGHS.Optimizer     # Solver
q         = 1 / S .* ones(S)    # Probabilidade dos cenarios
pu_Money  = 1e6;

# Leitura do PLD Newave
results_path = "March2024/Results/Solar_Eolica/"
dir = dirname(dirname(@__FILE__))
PLD = Matrix(DataFrame(CSV.File("March2024/PLD.csv";header=false)))'[:,:];

#Leitura número de horas no mês
h = Int.(Matrix(DataFrame(CSV.File("March2024/horas.csv";header=false)))[:,1]);


# ----------------------- CADASTRO USINAS -------------------------------

nI             = 2 # Número de usinas do caso
gu             = zeros(S, T, nI) # Geração das usinas
C              = zeros(nI,T)

GF_dis         = 50.0
GF_arinos      = 50.0
GF             = [GF_dis, GF_arinos]

TotGF = sum(GF);

#Leitura geração
gu_dis = Matrix(DataFrame(CSV.File("March2024/OUT-gen2.csv";header=true)))'[:,1:T];
gu_dis *=30
gu_dis /= mean(gu_dis)

gu_arinos = Matrix(DataFrame(CSV.File("March2024/OUT-gen1.csv";header=true)))'[:,1:T];
gu_arinos *=30;
gu_arinos /= mean(gu_arinos);

gu = ones(2000,T,nI);
gu[:,:,1] = gu_dis;
gu[:,:,2] = gu_arinos;

ymax = ones(nI); # Porcentagem máxima da geração das usinas contratadas
ymin = ones(nI); # Porcentagem minima da geração das usinas contratada


# ----------------------- CADASTRO CONTRATOS -------------------------------

J           = 1                         # Número de contratos
Qmax        = sum(GF[i] for i = 1:nI);   # Definição quantidade máxima
Q           = ones(J,T);                 # Vetor Q
Q[1,:]     .= Qmax;                      # Contrato A+1
if J > 1
    Q[2:13,:]  .= LinearAlgebra.I(12).*Qmax; # Contratos M+1
end

P           = ones(J, T);                # Vetor P

#a_Y         = 1.0;                      # Coeficiente angular p/ definição do preço
a_Y         = 0.93;
#b_Y         = 0.1; 
b_Y         = 11.31;                     # Coeficiente linear p/ definição do preço
P[1,:]     .= a_Y*mean(PLD[1:S,1:T]) + b_Y;

PLD_mensal  = Statistics.mean(PLD, dims = 1);
a_M         = a_Y;
b_M         = b_Y;
if J > 1
    P[2:end, :] .= (a_M .* PLD_mensal .+ b_M) .* LinearAlgebra.I(T)
end


global gammax = 0.5;
xmax     = ones(J)*gammax;      # Porcentagem maxima de contratação mensal
xmin     = ones(J).*-gammax;    # Porcentagem mínima de contratação mensal
xmax[1]  = 1.0;                 # Porcentagem maxima de contratação anual
xmin[1]  = 0.0;                 # Porcentagem minima de contratação anual


# ----------------------- PARÂMETROS DE RISCO -------------------------------

α                = 0.95 # Percentil CVaR Risco
flag_restricao   = false # Boolean restrição receita minima
Rmin             = ones(T)*(0.0); # Receita minima



# =========================================================================== #
#                               Run Omega Model                               #
# =========================================================================== #

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

# ----------------------- SOLVING PROPOSED MODEL -------------------------------

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



# ----------------------- MODELO ***UPSIDE*** TREINADO - Compute Statistics -------------------------------


Set_Quant       = [0.01 ; 0.05 ; 0.10 ; 0.25 ; 0.50 ; 0.75 ; 0.90 ; 0.95 ; 0.99];
nQuant          = length(Set_Quant);
Quant_Prop      = zeros(1, nQuant);

for (iter_q, q_) ∈ enumerate(Set_Quant)
    Quant_Prop[iter_q] = sort(R_otimo_proposto)[Int(round(q_*S))];
end

Avg_Prop        = mean(R_otimo_proposto);
CVaR_Prop       = Statistics.mean(sort(R_otimo_proposto)[1:Int(round((1 - α) * S)), :]);
CVaRUp_Prop       = Statistics.mean(sort(R_otimo_proposto)[Int(round((α) * S))+1:S,:]);


# ----------------------- MODELO ***UPSIDE*** TESTE - Compute Statistics -------------------------------

S_Test         = 800
R_teste_upside = zeros(S_Test)

for s in 1:S_Test
    R_teste_upside[s] = sum(sum(gu[s+S, t, i] * GF[i] * PLD[s+S, t] * h[t] for i in 1:nI) 
                        + sum((P[j, t] - PLD[s+S, t]) * Q[j, t] * x_novo[j] * h[t] for j in 1:J) for t in 1:T)./1e6;
end

Quant_Prop_Test      = zeros(1, nQuant);

for (iter_q, q_) ∈ enumerate(Set_Quant)
    Quant_Prop_Test[iter_q] = sort(R_teste_upside)[Int(round(q_*S_Test))];
end

Avg_Prop_Test        = mean(R_teste_upside);
CVaR_Prop_Test       = Statistics.mean(sort(R_teste_upside)[1:Int(round((1 - α) * S_Test)), :]);
CVaRUp_Prop_Test       = Statistics.mean(sort(R_teste_upside)[Int(round((α) * S_Test))+1:S_Test,:]);




# =========================================================================== #
#                               Run Risk-Reward Model                         #
# =========================================================================== #


Set_Λ           = [0.00 ; 0.05 ; 0.10 ; 0.25 ; 0.50 ; 0.75 ; 0.90 ; 0.95 ; 0.99];
n_Λ             = length(Set_Λ);

global Rdisp_RR        = zeros(n_Λ, S, T, nI);
global Rquant_RR       = zeros(n_Λ, S, T, J);
global RPort_RR        = zeros(n_Λ, S, T);
global RPortTot_RR     = zeros(n_Λ, S);
global xOptimal_RR     = zeros(n_Λ, J);
global yOptimal_RR     = zeros(n_Λ, nI);

Quant_RR        = zeros(n_Λ, nQuant);
Avg_RR          = zeros(n_Λ);
CVaR_RR         = zeros(n_Λ);
CVaRUp_RR       = zeros(n_Λ);

for (iter_λ, λ_) ∈ enumerate(Set_Λ)
    println("\n");
    println("λ = ", λ_);
    println("\n");

    global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
    Rdisp_RR[iter_λ, :, :, :], Rquant_RR[iter_λ, :, :, :], RPort_RR[iter_λ, :, :], RPortTot_RR[iter_λ, :], xOptimal_RR[iter_λ, :], yOptimal_RR[iter_λ, :] = port_alloc(params);
    
    # ----------------------- MODELO ***RISK-REWARD*** TREINADO - Compute Statistics -------------------------------

    for (iter_q, q_) ∈ enumerate(Set_Quant)
        Quant_RR[iter_λ, iter_q]    = sort(RPortTot_RR[iter_λ, :])[Int(round(q_*S))]/pu_Money;
        Avg_RR[iter_λ]              = mean(RPortTot_RR[iter_λ, :])/pu_Money;
        CVaR_RR[iter_λ]             = Statistics.mean(sort(RPortTot_RR[iter_λ, :])[1:Int(round((1 - α) * S)), :])/pu_Money;
        CVaRUp_RR[iter_λ]           = Statistics.mean(sort(RPortTot_RR[iter_λ, :])[Int(round((α) * S))+1:S, :])/pu_Money;

    end;

end;


# ----------------------- MODELO ***RISK-REWARD*** TESTE - Compute Statistics -------------------------------


S_Test  = 800
R_teste_RR = zeros(n_Λ, S_Test)

Quant_RR_Test        = zeros(n_Λ, nQuant);
Avg_RR_Test          = zeros(n_Λ);
CVaR_RR_Test         = zeros(n_Λ);
CVaRUp_RR_Test       = zeros(n_Λ);

for (iter_λ, λ_) ∈ enumerate(Set_Λ)
    println("\n");
    println("λ = ", λ_);
    println("\n");

    for s in 1:S_Test
        R_teste_RR[iter_λ, s] = sum(sum(gu[s+S, t, i] * GF[i] * PLD[s+S, t] * h[t] for i in 1:nI) 
                        + sum((P[j, t] - PLD[s+S, t]) * Q[j, t] * xOptimal_RR[iter_λ, j] * h[t] for j in 1:J) for t in 1:T);
    end

    for (iter_q, q_) ∈ enumerate(Set_Quant)
        Quant_RR_Test[iter_λ, iter_q]    = sort(R_teste_RR[iter_λ, :])[Int(round(q_*S_Test))]/pu_Money;
        Avg_RR_Test[iter_λ]              = mean(R_teste_RR[iter_λ, :])/pu_Money;
        CVaR_RR_Test[iter_λ]             = Statistics.mean(sort(R_teste_RR[iter_λ, :])[1:Int(round((1 - α) * S_Test)), :])/pu_Money;
        CVaRUp_RR_Test[iter_λ]           = Statistics.mean(sort(R_teste_RR[iter_λ, :])[Int(round((α) * S_Test))+1:S_Test, :])/pu_Money;

    end;

end;


# ================================================================= #

# ================================================================= #
#                      Structuring the Results                      #
# ================================================================= #

# ----------------------- Create DataFrames TREINADO -------------------------------

row_titles = ["0.00" ; "0.05" ; "0.10" ; "0.25" ; "0.50" ; "0.75" ; "0.90" ; "0.95" ; "0.99"; "Proposed"]
xMat    = round.([xOptimal_RR ; x_novo'].*100, digits = 2);
x_df    = DataFrame(xMat, :auto);
insertcols!(x_df, 1, "λ" => row_titles)

column_titles = ["Avg", "CVaR", "CVaR Up", "q 1%", "q 5%", "q 10%", "q 25%", "q 50%", "q 75%", "q 90%", "q 95%", "q 99%"]

println("\n\n");
println(" --> Contracting Levels <-- ");
println("\n");
print(x_df);


qMat    = [Quant_RR ; Quant_Prop];
AvgMat  = [Avg_RR ; Avg_Prop];
CVaRMat = [CVaR_RR ; CVaR_Prop];
CVaRUpMat = [CVaRUp_RR ; CVaRUp_Prop];

qMat    = [AvgMat CVaRMat CVaRUpMat qMat];
q_df    = DataFrame(qMat, :auto);
q_df = DataFrame(qMat, Symbol.(column_titles))

println("\n\n");
println(" --> Quantile Levels - in-sample <-- ");
println("\n");
print(q_df);
println("\n\n\n");


# ----------------------- Create DataFrames TESTE -------------------------------


qMat_Test    = [Quant_RR_Test ; Quant_Prop_Test];
AvgMat_Test  = [Avg_RR_Test ; Avg_Prop_Test];
CVaRMat_Test = [CVaR_RR_Test ; CVaR_Prop_Test];
CVaRUpMat_Test = [CVaRUp_RR_Test ; CVaRUp_Prop_Test];

qMat_Test    = [AvgMat_Test CVaRMat_Test CVaRUpMat_Test qMat_Test];
q_df_Test    = DataFrame(qMat_Test, :auto);

q_df_Test = DataFrame(qMat_Test, Symbol.(column_titles))

println("\n\n");
println(" --> Quantile Levels - out-of-sample <-- ");
println("\n");
print(q_df_Test);
println("\n\n\n");

# ----------------------- Indicador Risco-Upside -------------------------------

pos_Neut            = 1;
pos_Aver            = 5;
pos_Prop            = 10;
pos_CVaR            = 2;
pos_CVaRUp          = 3;

Indicador_InSample = (q_df[pos_Prop, pos_CVaRUp] - q_df[pos_Aver, pos_CVaRUp])/
                     (q_df[pos_Prop, pos_CVaR] - q_df[pos_Aver, pos_CVaR])

Indicador_OutofSample = (q_df_Test[pos_Prop, pos_CVaRUp] - q_df_Test[pos_Aver, pos_CVaRUp])/
                     (q_df_Test[pos_Prop, pos_CVaR] - q_df_Test[pos_Aver, pos_CVaR])

Indicators = round.([Indicador_InSample, Indicador_OutofSample], digits=2)
titles = ["In Sample", "Out of Sample"]
df_Indicator = DataFrame(Type = titles,Indicator = Indicators)

println("\n\n");
println(" --> Performance Indicator - In-sample <-- ");
println("\n");
print(df_Indicator);
println("\n\n\n");

# ------------------------ CSV Resultados -------------------------------------

xlsx_file_path = results_path*"Results_Stats.xlsx"

XLSX.openxlsx(xlsx_file_path, mode="w") do xf

    sheet_x_df = XLSX.addsheet!(xf, "x_df")
    XLSX.writetable!(sheet_x_df, x_df)

    sheet_q_df = XLSX.addsheet!(xf, "q_df")
    XLSX.writetable!(sheet_q_df, q_df)
    
    sheet_q_df_Test = XLSX.addsheet!(xf, "q_df_Test")
    XLSX.writetable!(sheet_q_df_Test, q_df_Test)

    sheet_df_Indicator = XLSX.addsheet!(xf, "Indicator")
    XLSX.writetable!(sheet_df_Indicator, df_Indicator)

end


# ================================================================= #

# ================================================================= #
#                      Plots                                        #
# ================================================================= #

#
    # ---> Plot: Inverse Cumulative Distribution <---
#

s_linewidth         = 2;
λ_plot              = Set_Λ[pos_Aver];
delta               = 10;
α_d                 = 0.20;
name                = "InSample"

inv_cum_graph(s_linewidth, pos_Neut, pos_Aver, λ_plot, delta,
                R_otimo_proposto, RPortTot_RR, pu_Money, α_d, S, name)

name = "OutOfSample"
inv_cum_graph(s_linewidth, pos_Neut, pos_Aver, λ_plot, delta,
                R_teste_upside, R_teste_RR, pu_Money, α_d, S_Test, name)


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

contrato[:,1]   .= soma_neutral*Qmax./100;
contrato[:,2]   .= soma_averse*Qmax./100;
contrato[:,3]   .= soma_proposto*Qmax./100;

# -------------- InSample Plot ----------------

geracao1_InSample        = gu_dis[1:S,:]*GF[1];
geracao2_InSample        = gu_arinos[1:S,:]*GF[2];

p_energy_InSample = plot();

p_energy_InSample = groupedbar!(
    contrato,
    color = ["black" "gray" "blue"],
    linecolor = :white,
    fillalpha = [0.9 0.9 0.9],
    label = ["Risk-Neutral" string("Risk-Averse (λ = ", λ_plot, ")") "Proposed Approach"],
    xlabel = "Bussiness Periods (Months)",
    ylabel = "Energy (avgMW)",
    ylims = (50,300),
    legend = :bottomright,
    linewidth = 1,
    x_foreground_color_border = :white
);

xaxis_leg       = [i for i ∈ 1:12];
vRES_Generation = zeros(S,T)
for s ∈ 1:S
    for t ∈ 1:T
        vRES_Generation[s,t] = (geracao1_InSample[s,t] + geracao2_InSample[s,t]);
    end;
end;

lev_fillalpha   = 0.5;
lev_color       = :limegreen;

p_energy_InSample.n = 0;

p_energy_InSample = violin!(
    vRES_Generation,
    color = lev_color,
    ylims = (0,180),
    legend = :topleft,
    xticks = 1:1:12,
    x_foreground_color_border = :white,
    xlabel = "Bussiness Periods (Months)",
    fillalpha=[lev_fillalpha lev_fillalpha lev_fillalpha],
    label = ["Renewable Production" false false false false false false false false false false false false],
    linecolor = nothing
);

display(p_energy_InSample);
Plots.savefig(p_energy_InSample, results_path*"p_energy_InSample.png");

# -------------- OutOfSample Plot ----------------

geracao1_OutOfSample        = gu_arinos[S+1:S+S_Test,:]*GF[1];
geracao2_OutOfSample        = gu_dis[S+1:S+S_Test,:]*GF[2];

p_energy_OutOfSample = plot();

p_energy_OutOfSample = groupedbar!(
    contrato,
    color = ["black" "gray" "blue"],
    linecolor = :white,
    fillalpha = [0.9 0.9 0.9],
    label = ["Risk-Neutral" string("Risk-Averse (λ = ", λ_plot, ")") "Proposed Approach"],
    xlabel = "Bussiness Periods (Months)",
    ylabel = "Energy (avgMW)",
    ylims = (50,300),
    legend = :bottomright,
    linewidth = 1,
    x_foreground_color_border = :white
);

xaxis_leg       = [i for i ∈ 1:12];
vRES_Generation = zeros(S_Test,T)
for s ∈ 1:S_Test
    for t ∈ 1:T
        vRES_Generation[s,t] = (geracao1_OutOfSample[s,t] + geracao2_OutOfSample[s,t]);
    end;
end;

lev_fillalpha   = 0.5;
lev_color       = :limegreen;

p_energy_OutOfSample.n = 0;

p_energy_OutOfSample = violin!(
    vRES_Generation,
    color = lev_color,
    ylims = (0,180),
    legend = :topleft,
    xticks = 1:1:12,
    x_foreground_color_border = :white,
    xlabel = "Bussiness Periods (Months)",
    fillalpha=[lev_fillalpha lev_fillalpha lev_fillalpha],
    label = ["Renewable Production" false false false false false false false false false false false false],
    linecolor = nothing
);

display(p_energy_OutOfSample);
Plots.savefig(p_energy_OutOfSample, results_path*"p_energy_OutOfSample.png");