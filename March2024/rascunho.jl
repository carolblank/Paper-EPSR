using PortfolioOptCTG
using Plots
using JuMP
using XLSX, DataFrames, Statistics, CSV
using StatsPlots
using LinearAlgebra
using HiGHS


param_caso = XLSX.readxlsx("Monthly_Report/August2023/Parametros.xlsx")["Parâmetros Caso"]

S = param_caso["F3"]; # Numero de cenários
T = param_caso["F6"]; # Numero de meses
tx = 0; # Custo oportunidade anual
tx_m = (1+tx)^(1/12) -1; # Custo oportunidade mensal
Solv = HiGHS.Optimizer # Solver
q = 1 / S .* ones(S); # Probabilidade dos cenarios

# Leitura do PLD Newave
file_PLD = string("Monthly_Report/August2023//",param_caso["E9"])
PLD = Matrix(DataFrame(CSV.File(file_PLD;header=false)))'[:,:];

#Leitura número de horas no mês
file_h = string("Monthly_Report/August2023//",param_caso["E12"])
h = Int.(Matrix(DataFrame(CSV.File(file_h;header=false)))[:,1]);

param_usinas = XLSX.readxlsx("Monthly_Report/August2023/Parametros.xlsx")["Cadastro Usinas"]

usina = Vector() # Usinas do caso

if param_usinas["D5"] == 1
      push!(usina, "ARINOS")
end

if param_usinas["E5"] == 1
      push!(usina, "DIS")
end

if param_usinas["F5"] == 1
      push!(usina, "SDP")
end

if param_usinas["G5"] == 1
      push!(usina, "HYDRO")
end

nI = length(usina) # Número de usinas do caso
gu = zeros(S, T, nI) # Geração das usinas
C = ones(nI, T) # Custo de construção/manutenção por usina por mês
n = 1
GF = ones(nI)
ymax = ones(nI); # Porcentagem máxima da geração das usinas contratadas
ymin = zeros(nI); # Porcentagem minima da geração das usinas contratada


#Leitura geração usina ARINOS (% garantia física)
arinos = string("Monthly_Report/August2023/",param_usinas["D11"])
gu_arinos = Matrix(DataFrame(CSV.File(arinos;header=false)))'[1:S,1:T];

#Leitura geração usina DIS (% garantia física)
dis = string("Monthly_Report/August2023//",param_usinas["E11"])
gu_dis = Matrix(DataFrame(CSV.File(dis;header=false)))'[1:S,1:T];

#Leitura geração usina SDP (% garantia física)
sdp = string("Monthly_Report/August2023/",param_usinas["F11"])
gu_sdp = Matrix(DataFrame(CSV.File(sdp;header=false)))'[1:S,1:T];

#Leitura geração usina DIS (% garantia física)
hydro = string("Monthly_Report/August2023/",param_usinas["G11"])
gu_hydro = Matrix(DataFrame(CSV.File(hydro;header=false)))'[1:S,1:T];

if param_usinas["D5"] == 1
      C[n,:] .= param_usinas["D9"]
      GF[n] = param_usinas["D17"]
      gu[:,:,n] .= gu_arinos
      ymax[n] = param_usinas["D21"]
      ymin[n] = param_usinas["D25"]
      n +=1
end

if param_usinas["E5"] == 1
      C[n,:] .= param_usinas["E9"]
      GF[n] = param_usinas["E17"]
      gu[:,:,n] .= gu_dis
      ymax[n] = param_usinas["E21"]
      ymin[n] = param_usinas["E25"]
      n +=1
end

if param_usinas["F5"] == 1
      C[n,:] .= param_usinas["F9"]
      GF[n] = param_usinas["F17"]
      gu[:,:,n] .= gu_sdp
      ymax[n] = param_usinas["F21"]
      ymin[n] = param_usinas["F25"]
      n +=1
end

if param_usinas["G5"] == 1
      C[n,:] .= param_usinas["G9"]
      GF[n] = param_usinas["G17"]
      gu[:,:,n] .= gu_hydro
      ymax[n] = param_usinas["G21"]
      ymin[n] = param_usinas["G25"]
      n +=1
end

param_contrato = XLSX.readxlsx("Monthly_Report/August2023/Parametros.xlsx")["Cadastro Contrato Quantidade"]

Q = ones(0,T)
P = ones(0,T)
xmax = ones(0)
xmin = zeros(0)
n=1

if param_contrato["D5"] == 1       # Contrato 1
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F8:Q8"]
     P[n,:] = param_contrato["F11:Q11"]
     xmin[n] = param_contrato["T7"]
     xmax[n] = param_contrato["T10"]
     n+=1
end

if param_contrato["D16"] == 1      # Contrato 2
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F19:Q19"]
     P[n,:] = param_contrato["F22:Q22"]
     xmin[n] = param_contrato["T18"]
     xmax[n] = param_contrato["T21"]
     n+=1
end

if param_contrato["D24"] == 1      # Contrato 3
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F27:Q27"]
     P[n,:] = param_contrato["F30:Q30"]
     xmin[n] = param_contrato["T26"]
     xmax[n] = param_contrato["T29"]
     n+=1
end

if param_contrato["D32"] == 1      # Contrato 4
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F35:Q35"]
     P[n,:] = param_contrato["F38:Q38"]
     xmin[n] = param_contrato["T34"]
     xmax[n] = param_contrato["T37"]
     n+=1
end

if param_contrato["D40"] == 1      # Contrato 5
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F43:Q43"]
     P[n,:] = param_contrato["F46:Q46"]
     xmin[n] = param_contrato["T42"]
     xmax[n] = param_contrato["T45"]
     n+=1
end

if param_contrato["D48"] == 1      # Contrato 6
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F51:Q51"]
     P[n,:] = param_contrato["F54:Q54"]
     xmin[n] = param_contrato["T50"]
     xmax[n] = param_contrato["T53"]
     n+=1
end

if param_contrato["D56"] == 1      # Contrato 7
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F59:Q59"]
     P[n,:] = param_contrato["F62:Q62"]
     xmin[n] = param_contrato["T58"]
     xmax[n] = param_contrato["T61"]
     n+=1
end

if param_contrato["D64"] == 1      # Contrato 8
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F67:Q67"]
     P[n,:] = param_contrato["F70:Q70"]
     xmin[n] = param_contrato["T66"]
     xmax[n] = param_contrato["T69"]
     n+=1
end

if param_contrato["D72"] == 1      # Contrato 9
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F75:Q75"]
     P[n,:] = param_contrato["F78:Q78"]
     xmin[n] = param_contrato["T74"]
     xmax[n] = param_contrato["T77"]
     n+=1
end

if param_contrato["D80"] == 1      # Contrato 10
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F83:Q83"]
     P[n,:] = param_contrato["F86:Q86"]
     xmin[n] = param_contrato["T82"]
     xmax[n] = param_contrato["T85"]
     n+=1
end

if param_contrato["D88"] == 1      # Contrato 11
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F91:Q91"]
     P[n,:] = param_contrato["F94:Q94"]
     xmin[n] = param_contrato["T90"]
     xmax[n] = param_contrato["T93"]
     n+=1
end

if param_contrato["D96"] == 1      # Contrato 12
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F99:Q99"]
     P[n,:] = param_contrato["F102:Q102"]
     xmin[n] = param_contrato["T98"]
     xmax[n] = param_contrato["T101"]
     n+=1
end

if param_contrato["D104"] == 1      # Contrato 13
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F107:Q107"]
     P[n,:] = param_contrato["F110:Q110"]
     xmin[n] = param_contrato["T106"]
     xmax[n] = param_contrato["T109"]
     n+=1
end

if param_contrato["D118"] == 1      # Contrato 14
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F121:Q121"]
     P[n,:] = param_contrato["F124:Q124"]
     xmin[n] = param_contrato["T120"]
     xmax[n] = param_contrato["T123"]
     n+=1
end

if param_contrato["D118"] == 1      # Contrato 15
     Q = vcat(Q, ones(T)')
     P = vcat(P, ones(T)')
     xmax = vcat(xmax,1)
     xmin = vcat(xmin,1)
     Q[n,:] = param_contrato["F132:Q132"]
     P[n,:] = param_contrato["F135:Q135"]
     xmin[n] = param_contrato["T131"]
     xmax[n] = param_contrato["T134"]
     n+=1
end

J = n-1

println(Q)
println(P)
println(xmax)
println(xmin)

param_risco = XLSX.readxlsx("Monthly_Report/August2023/Parametros.xlsx")["Parâmetros de Risco"]

α = 0.01*param_risco["F3"] # Percentil CVaR Risco
flag_restricao = param_risco["F6"] # Boolean restrição receita minima
λ = 0.01*param_risco["F9"] # Reguçador de risco função objetivo
Rmin = ones(T)*(1e6)*param_risco["F12"]; # Receita minima

params = ParametrosContratacao(; gu, GF, h, Q, PLD, α, T, λ, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv)
Rdisp_otimo, Rquant_otimo, R_otimo, R_total_otimo, x_otimo, y_otimo = port_alloc(params);


x = zeros(15)
y = zeros(4)

m=1
if param_usinas["D5"] == 1
    println("Construção ótima ", usina[m],": ",GF[m]*y_otimo[m],"MW méd")
    y[1] = y_otimo[m]
    m+=1
else
    y[1] = 0
end
if param_usinas["E5"] == 1
    println("Construção ótima ", usina[2],": ",GF[2]*y_otimo[2],"MW méd")
    y[2] = y_otimo[m]
    m+=1
else
    y[2] = 0
end
if param_usinas["F5"] == 1
    println("Construção ótima ", usina[3],": ",GF[3]*y_otimo[3],"MW méd")
    y[3] = y_otimo[m]
    m+=1
else
    y[3] = 0
end
if param_usinas["G5"] == 1
    println("Construção ótima ", usina[4],": ",GF[4]*y_otimo[4],"MW méd")
    y[4] = y_otimo[m]
else
    y[4] = 0
end

n=1
if param_contrato["D5"] == 1
    println("Quantidade ótima A+1: ", Q[n]*x_otimo[n],"MW méd")
    x[1] = x_otimo[n]
    n+=1
else
    x[1] = 0
end
if param_contrato["D16"] == 1 
    println("Quantidade ótima M+1 - Jan: ", Q[n,1]*x_otimo[n],"MW méd")
    x[2] = x_otimo[n]
    n+=1
else
    x[2] = 0
end
if param_contrato["D24"] == 1 
    println("Quantidade ótima M+1 - Fev: ", Q[n,2]*x_otimo[n],"MW méd")
    x[3] = x_otimo[n]
    n+=1
else
    x[3] = 0
end
if param_contrato["D32"] == 1
    println("Quantidade ótima M+1 - Mar: ", Q[n,3]*x_otimo[n],"MW méd")
    x[4] = x_otimo[n]
    n+=1
else
    x[4] = 0
end
if param_contrato["D40"] == 1
    println("Quantidade ótima M+1 - Abr: ", Q[n,4]*x_otimo[n],"MW méd")
    x[5] = x_otimo[n]
    n+=1
else
    x[5] = 0
end
if param_contrato["D48"] == 1
    println("Quantidade ótima M+1 - Mai: ", Q[n,5]*x_otimo[n],"MW méd")
    x[6] = x_otimo[n]
    n+=1
else
    x[6] = 0
end
if param_contrato["D56"] == 1
    println("Quantidade ótima M+1 - Jun: ", Q[n,6]*x_otimo[n],"MW méd")
    x[7] = x_otimo[n]
    n+=1
else
    x[7] = 0
end
if param_contrato["D64"] == 1
    println("Quantidade ótima M+1 - Jul: ", Q[n,7]*x_otimo[n],"MW méd")
    x[8] = x_otimo[n]
    n+=1
else
    x[8] = 0
end
if param_contrato["D72"] == 1
    println("Quantidade ótima M+1 - Ago: ", Q[n,8]*x_otimo[n],"MW méd")
    x[9] = x_otimo[n]
    n+=1
else
    x[9] = 0
end
if param_contrato["D80"] == 1
    println("Quantidade ótima M+1 - Set: ", Q[n,9]*x_otimo[n],"MW méd")
    x[10] = x_otimo[n]
    n+=1
else
    x[10] = 0
end
if param_contrato["D88"] == 1
    println("Quantidade ótima M+1 - Out: ", Q[n,10]*x_otimo[n],"MW méd")
    x[11] = x_otimo[n]
    n+=1
else
    x[11] = 0
end
if param_contrato["D96"] == 1
    println("Quantidade ótima M+1 - Nov: ", Q[n,11]*x_otimo[n],"MW méd")
    x[12] = x_otimo[n]
    n+=1
else
    x[12] = 0
end
if param_contrato["D104"] == 1
    println("Quantidade ótima M+1 - Dez: ", Q[n,12]*x_otimo[n],"MW méd")
    x[13] = x_otimo[n]
    n+=1
else
    x[13] = 0
end
if param_contrato["D118"] == 1
    println("Quantidade ótima - Genérico: ", Q[n]*x_otimo[n],"MW méd")
    x[14] = x_otimo[n]
    n+=1
else
    x[14] = 0
end
if param_contrato["D129"] == 1
    println("Quantidade ótima - Genérico: ", Q[n]*x_otimo[n],"MW méd")
    x[15] = x_otimo[n]
else
    x[15] = 0
end

println("% Contratação: ", x_otimo)

out_Receita = DataFrame(R_otimo',:auto);

XLSX.openxlsx("Monthly_Report/August2023/Parametros.xlsx", mode="rw") do xf
    sheet = xf[5]
    sheet["A2"] = Statistics.mean(R_otimo'/1e6, dims = 2)
end

x_otimo