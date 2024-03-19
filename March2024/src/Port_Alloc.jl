using JuMP
using XLSX, DataFrames, Statistics
using CSV
using LinearAlgebra
using HiGHS

"""
    ParametrosContratacao{Z<:Real}

Parametros para a alocacao de portfolio.

Parametros:
    - gu::Array{Z, 3} # Geração das usinas
    - GF::Vector{Z} # Garatia fisica
    - h::Vector{Int} # numero de horas 
    - Q::Array{Z, 2} # Quantidade maxima para contratação
    - PLD::Array{Z, 2} # Preço
    - α::Z # Percentil CVaR Risco
    - T::Int # Numero de meses
    - λ::Z # Reguçador de risco função objetivo
    - flag_restricao::Bool # Boolean restrição receita minima
    - Rmin::Vector{Z} # Receita minima
    - P::Array{Z, 2} # Preço contrato futuro
    - S::Int # numero de cenários 
    - q::Vector{Z} # probabilidade dos cenarios
    - nI::Int # Número de usinas disponíveis
    - J::Int # Número de contratos
    - C::Array{Z,2} # Custo de construção/manutenção por usina por mês
    - tx_m::Z # Custo oportunidade mensal
    - xmin::Vector{Z} # porcentagem minima de contratação
    - xmax::Vector{Z} # porcentagem maxima de contratação 
    - ymin::Vector{Z} # porcentagem minima da geração das usinas contratada
    - ymax::Vector{Z} # porcentagem minima da geração das usinas contratada
    - Solv # solver
"""
struct ParametrosContratacao{Z<:Real}
    gu::Array{Z,3} # Geração das usinas
    GF::Vector{Z} # Garatia fisica
    h::Vector{Int} # numero de horas 
    Q::Array{Z,2} # Quantidade maxima para contratação
    PLD::Array{Z,2} # Preço
    α::Z # Quantil CVar Risco
    T::Int # Numero de meses
    λ::Z # Reguçador de risco função objetivo
    flag_restricao::Bool # Boolean restrição receita minima
    Rmin::Union{Vector{Z},Nothing} # Receita minima
    P::Array{Z,2} # Preço contrato futuro
    S::Int # numero de cenários 
    q::Vector{Z} # probabilidade dos cenarios
    nI::Int # Número de usinas disponíveis
    J::Int # Número de contratos
    C::Array{Z,2} # Custo de construção/manutenção por usina por mês
    tx_m::Z # Custo oportunidade mensal
    xmin::Vector{Z} # porcentagem minima de contratação
    xmax::Vector{Z} # porcentagem maxima de contratação 
    ymin::Vector{Z} # porcentagem minima da geração das usinas contratada
    ymax::Vector{Z} # porcentagem minima da geração das usinas contratada
    Solv::Any # solver

    function ParametrosContratacao(
        gu::Array{Z,3},
        GF::Vector{Z},
        h::Vector{Int},
        Q::Array{Z,2},
        PLD::Array{Z,2},
        α::Z,
        T::Int,
        λ::Z,
        flag_restricao::Bool,
        Rmin::Union{Vector{Z},Nothing},
        P::Array{Z,2},
        S::Int,
        q::Vector{Z},
        nI::Int,
        J::Int,
        C::Array{Z,2},
        tx_m::Z,
        xmin::Vector{Z},
        xmax::Vector{Z},
        ymin::Vector{Z},
        ymax::Vector{Z},
        Solv::Any
    ) where {Z<:Real}
        α >= 0 || throw(ArgumentError("α must be >= 0"))
        α <= 1 || throw(ArgumentError("α must be <= 1"))

        return new{Z}(
            gu,
            GF,
            h,
            Q,
            PLD,
            α,
            T,
            λ,
            flag_restricao,
            Rmin,
            P,
            S,
            q,
            nI,
            J,
            C,
            tx_m,
            xmin,
            xmax,
            ymin,
            ymax,
            Solv
        )
    end
end

function ParametrosContratacao(; gu, GF, h, Q, PLD, α, T, λ, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv)

    return ParametrosContratacao(
        gu,
        GF,
        h,
        Q,
        PLD,
        α,
        T,
        λ,
        flag_restricao,
        Rmin,
        P,
        S,
        q,
        nI,
        J,
        C,
        tx_m,
        xmin,
        xmax,
        ymin,
        ymax,
        Solv
    )
    
end

"""
    port_alloc(params::ParametrosContratacao)

Funcao para Allocacao de portifolio.
"""
function port_alloc(params::ParametrosContratacao)

    gu,
    GF,
    h,
    Q,
    PLD,
    α,
    T,
    λ,
    flag_restricao,
    Rmin,
    P,
    Solv,
    S,
    q,
    nI,
    J,
    C,
    tx_m,
    xmin,
    xmax,
    ymin,
    ymax = params.gu,
    params.GF,
    params.h,
    params.Q,
    params.PLD,
    params.α,
    params.T,
    params.λ,
    params.flag_restricao,
    params.Rmin,
    params.P,
    params.Solv,
    params.S,
    params.q,
    params.nI,
    params.J,
    params.C,
    params.tx_m,
    params.xmin,
    params.xmax,
    params.ymin,
    params.ymax

    model = Model(Solv)

    @variable(model, R_disp[1:S, 1:T, 1:nI])
    @variable(model, R_quant[1:S, 1:T, 1:J])
    @variable(model, R[1:S, 1:T])
    @variable(model, w1)
    @variable(model, δ1[1:S] >= 0)
    @variable(model, x[1:J])
    @variable(model, y[1:nI])
    @variable(model, R_total[1:S])

    @objective(
        model,
        Max,
        λ * (w1 - sum(δ1[s] * q[s] / (1 - α) for s = 1:S)) +
        (1 - λ) * (sum(q[s] * R_total[s] for s = 1:S))
    )

    @constraint(model, [s ∈ 1:S], R_total[s] == sum(R[s, t]/((1+tx_m)^t) for t = 1:T))

    @constraint(
        model,
        [s ∈ 1:S, t ∈ 1:T, i ∈ 1:nI],
        R_disp[s, t, i] ==
        gu[s, t, i] * GF[i] * y[i] * PLD[s, t] * h[t] - (C[i, t] * (GF[i] * y[i]) * h[t])
    )

    @constraint(
        model,
        [s ∈ 1:S, t ∈ 1:T, j ∈ 1:J],
        R_quant[s, t, j] == (P[j, t] - PLD[s, t]) * Q[j, t] * x[j] * h[t]
    )

    @constraint(
        model,
        [s ∈ 1:S, t ∈ 1:T],
        R[s, t] == sum(R_disp[s, t, i] for i = 1:nI) + sum(R_quant[s, t, j] for j = 1:J)
    )

    @constraint(
        model,
        [t ∈ 1:T],
        sum(Q[j, t] * x[j] for j = 1:J) <= sum(GF[i] * (1 + gammax) * y[i] for i = 1:nI)
    )

    @constraint(model, [j ∈ 1:J], xmin[j] <= x[j])

    @constraint(model, [j ∈ 1:J], x[j] <= xmax[j])

    @constraint(model, [i ∈ 1:nI], ymin[i] <= y[i])

    @constraint(model, [i ∈ 1:nI], y[i] <= ymax[i])

    @constraint(model, [s ∈ 1:S], δ1[s] >= w1 - R_total[s])

    if flag_restricao
        @variable(model, w2[1:T])
        @variable(model, δ2[1:S, 1:T] >= 0)

        @constraint(
            model,
            [t ∈ 1:T],
            w2[t] - sum(δ2[s, t] * q[s] / (1 - α) for s = 1:S) >= Rmin[t]
        )

        @constraint(model, [s ∈ 1:S, t ∈ 1:T], δ2[s, t] >= w2[t] - R[s, t])

    end

    status = optimize!(model)
    Rdisp_otimo = JuMP.value.(R_disp)
    Rquant_otimo = JuMP.value.(R_quant)
    x_otimo = JuMP.value.(x)
    y_otimo = JuMP.value.(y)
    R_otimo = JuMP.value.(R)
    R_total_otimo = JuMP.value.(R_total)
    resp = JuMP.objective_value(model)

    return Rdisp_otimo, Rquant_otimo, R_otimo, R_total_otimo, x_otimo, y_otimo
end


function inv_cum_graph(s_linewidth, pos_Neut, pos_Aver, λ_plot, delta,
                        R_upside, R_RR, pu_Money, α_d, S, name)

    Sd = Int(round(S*α_d));

    min_xlim = Int(round(min(minimum(R_upside), 
    minimum(R_RR[pos_Aver, :]./pu_Money), 
    minimum(R_RR[pos_Neut, :]./pu_Money)
    )/10, digits = 0)*10) - delta;

    max_xlim = Int(round(max(maximum(R_upside), 
    maximum(R_RR[pos_Aver, :]./pu_Money), 
    maximum(R_RR[pos_Neut, :]./pu_Money)
    )/10, digits = 0)*10) + delta;

    p1 = plot();

    R_Prop_Sort         = sort(R_upside);

    p1 = plot!((1:S)./S, R_Prop_Sort, ylims = (min_xlim,max_xlim), labels = "Proposed approach",
    legend = :topleft, color=:blue, xlabel="Inverse Cumulative Probability",
    xticks = 0.0:0.1:1.0,
    yticks = min_xlim:10:max_xlim,
    ylabel = "Revenue MMR\$",
    linewidth = s_linewidth
    );

    R_RR_Sort_Neut      = sort(R_RR[pos_Neut, :]);

    p1 = plot!((1:S)./S, R_RR_Sort_Neut./pu_Money, 
    color = :black, 
    labels = "Risk-Neutral", 
    linewidth = s_linewidth,
    line=(:dot, 3)
    );

    R_RR_Sort_Aver      = sort(R_RR[pos_Aver, :]);

    p1 = plot!((1:S)./S, R_RR_Sort_Aver./pu_Money, 
    color = :gray, 
    labels = string("Risk-Averse (λ = ", λ_plot, ")"),
    linewidth = s_linewidth,
    ls=:dash
    );

    display(p1);
    fig_title1 = string("InvDistAcum_Full_",name,".png")
    Plots.savefig(p1, fig_title1);

    #
    # ---> Plot: Inverse Cumulative Distribution -- Downside Region <---
    #

    delta = 5;

    min_xlim = Int(round(min(minimum(R_Prop_Sort[1:Sd]), 
    minimum(R_RR_Sort_Aver[1:Sd]./pu_Money), 
    minimum(R_RR_Sort_Neut[1:Sd]./pu_Money)
    )/10, digits = 0)*10) - delta;

    max_xlim = Int(round(max(maximum(R_Prop_Sort[1:Sd]), 
    maximum(R_RR_Sort_Aver[1:Sd]./pu_Money), 
    maximum(R_RR_Sort_Neut[1:Sd]./pu_Money)
    )/10, digits = 0)*10) + delta;

    p2 = plot();

    R_Prop_Sort         = sort(R_upside);

    x_axis__            = 0:0.2/(Sd-1):α_d;

    p2 = plot!(x_axis__, R_Prop_Sort[1:Sd], ylims = (min_xlim,max_xlim), labels = "Proposed approach",
    legend = :bottomright, color=:blue, xlabel="Inverse Cumulative Probability",
    xlims = (0.00,(α_d+0.01)),
    xticks = 0.00:0.02:α_d,
    ylabel = "Revenue MMR\$",
    linewidth = s_linewidth
    );

    R_RR_Sort_Neut      = sort(R_RR[pos_Neut, :]);

    p2 = plot!(x_axis__, R_RR_Sort_Neut[1:Sd]./pu_Money, 
    color = :black, 
    labels = "Risk-Neutral", 
    linewidth = s_linewidth,
    line=(:dot, 3)
    );

    R_RR_Sort_Aver      = sort(R_RR[pos_Aver, :]);

    p2 = plot!(x_axis__, R_RR_Sort_Aver[1:Sd]./pu_Money, 
    color = :gray, 
    labels = string("Risk-Averse (λ = ", λ_plot, ")"),
    linewidth = s_linewidth,
    ls=:dash
    );

    fig_title2 = string("InvDistAcum_DownSide_",name,".png")
    Plots.savefig(p2, fig_title2);

    #
    # ---> Plot: Inverse Cumulative Distribution -- Upside Region <---
    #

    p3 = plot();

    delta = 5;

    min_xlim = Int(round(min(minimum(R_Prop_Sort[S-Sd+1:S]), 
    minimum(R_RR_Sort_Aver[S-Sd+1:S]./pu_Money), 
    minimum(R_RR_Sort_Neut[S-Sd+1:S]./pu_Money)
    )/10, digits = 0)*10) - delta;

    max_xlim = Int(round(max(maximum(R_Prop_Sort[S-Sd+1:S]), 
    maximum(R_RR_Sort_Aver[S-Sd+1:S]./pu_Money), 
    maximum(R_RR_Sort_Neut[S-Sd+1:S]./pu_Money)
    )/10, digits = 0)*10) + delta;

    R_Prop_Sort         = sort(R_upside);
    x_axis__            = (1 - α_d):0.2/(Sd-1):1;


    p3 = plot!(x_axis__, R_Prop_Sort[S-Sd+1:S], 
    ylims = (min_xlim,max_xlim), 
    labels = "Proposed approach",
    legend = :topleft, color=:blue, xlabel="Inverse Cumulative Probability",
    xlims = ((1-α_d),1),
    xticks = (1-α_d):0.02:1,
    yticks = min_xlim:10:max_xlim,
    ylabel = "Revenue MMR\$",
    linewidth = s_linewidth
    );

    R_RR_Sort_Neut      = sort(R_RR[pos_Neut, :]);

    p3 = plot!(x_axis__, R_RR_Sort_Neut[S-Sd+1:S]./pu_Money, 
    color = :black, 
    labels = "Risk-Neutral", 
    linewidth = s_linewidth,
    line=(:dot, 3)
    );

    R_RR_Sort_Aver      = sort(R_RR[pos_Aver, :]);

    p3 = plot!(x_axis__, R_RR_Sort_Aver[S-Sd+1:S]./pu_Money, 
    color = :gray, 
    labels = string("Risk-Averse (λ = ", λ_plot, ")"),
    linewidth = s_linewidth,
    ls=:dash
    );

    fig_title3 = string("InvDistAcum_Upside_",name,".png")
    Plots.savefig(p3, fig_title3);

return

end