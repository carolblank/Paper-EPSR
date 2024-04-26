
# ----------------------------- Upside -------------------------------

range_P = collect(0:10:200)

x_otimo_disp_upside = zeros(length(range_P))
CVaR_Prop = zeros(length(range_P))
CVaRUp_Prop = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel

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

        x_otimo_disp_upside[iter] = x_novo[1]

        CVaR_Prop[iter]       = Statistics.mean(sort(R_otimo_proposto)[1:Int(round((1 - α) * S)), :]);
        CVaRUp_Prop[iter]      = Statistics.mean(sort(R_otimo_proposto)[Int(round((α) * S))+1:S,:]);

end

Y = x_otimo_disp_upside
X = range_P

p = plot(X, Y, 
        label = "Proposed approach",
        xlabel = "Contract Price (A+1)",
        ylabel = "Contracting Level",
        color=:blue,
        linewidth = s_linewidth)




# ----------------------------- Risk-Reward -------------------------------


x_otimo_disp_RR = zeros(length(range_P))
CVaR_RR = zeros(length(range_P))
CVaRUp_RR = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel

        λ_ = 0.5

        global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
        Rdisp_RR, Rquant_RR, RPort_RR, RPortTot_RR, xOptimal_RR, yOptimal_RR = port_alloc(params);
    
        x_otimo_disp_RR[iter] = xOptimal_RR[1]
        CVaR_RR[iter]             = Statistics.mean(sort(RPortTot_RR)[1:Int(round((1 - α) * S)), :])/pu_Money;
        CVaRUp_RR[iter]           = Statistics.mean(sort(RPortTot_RR)[Int(round((α) * S))+1:S, :])/pu_Money;

end

X = range_P
Y = x_otimo_disp_RR

p = plot!(X, Y, 
        label = "Risk-Averse (λ = 0.5)",
        color = :gray,
        linewidth = s_linewidth,
        ls=:dash)


# ----------------------------- Risk-neutral -------------------------------

x_otimo_disp_Neutral = zeros(length(range_P))
CVaR_Neutral = zeros(length(range_P))
CVaRUp_Neutral = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel

        λ_ = 0.0

        global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
        Rdisp_Neutral, Rquant_Neutral, RPort_Neutral, RPortTot_Neutral, xOptimal_Neutral, yOptimal_Neutral = port_alloc(params);
    
        x_otimo_disp_Neutral[iter] = xOptimal_Neutral[1]
        CVaR_Neutral[iter]             = Statistics.mean(sort(RPortTot_Neutral)[1:Int(round((1 - α) * S)), :])/pu_Money;
        CVaRUp_Neutral[iter]           = Statistics.mean(sort(RPortTot_Neutral)[Int(round((α) * S))+1:S, :])/pu_Money;

end


X = range_P
Y = x_otimo_disp_Neutral

p = plot!(X,Y, 
        label = "Risk-Neutral",
        color = :black,
        linewidth = s_linewidth,
        line=(:dot, 3))

display(p)

Plots.savefig(p, results_path*"Disp_Contratar.png");

plot_cvau = plot(X, CVaRUp_Prop - CVaRUp_RR, 
        label = "Risk-Averse",
        xlabel = "Contract Price (A+1)",
        ylabel = "CVAU",
        color = :gray,
        linewidth = s_linewidth,
        ls=:dash
);
plot_cvau = plot!(X, CVaRUp_Prop, 
        label = "Proposed approach",
        xlabel = "Contract Price (A+1)",
        ylabel = "CVAU",
        color = :blue,
        linewidth = s_linewidth
);
plot_cvau = plot!(X, CVaRUp_Prop - CVaRUp_Neutral, 
        label = "Risk-Neutral",
        color = :black,
        linewidth = s_linewidth,
        line=(:dot, 3)
)


Eta_RR = zeros(21)
CVaRUp_Compared_RR = zeros(21)
for i in 1:21
    Eta_RR[i] = (CVaRUp_Prop[i] - CVaRUp_RR[i])/(CVaR_Prop[i] - CVaR_RR[i])
end
for i in 1:21
    CVaRUp_Compared_RR[i] = (CVaRUp_Prop[i] - CVaRUp_RR[i])/CVaRUp_RR[i]
end
plot_eta = plot(X, Eta_RR, 
        label = "Risk-Averse (λ = 0.5)",
        xlabel = "Contract Price (A+1)",
        ylabel = "η (Perfomance Indicator)",
        color = :gray,
        linewidth = s_linewidth,
        ls=:dash
)

Eta_Neutral = zeros(21)
CVaRUp_Compared_Neutral = zeros(21)
for i in 1:21
    Eta_Neutral[i] = (CVaRUp_Prop[i] - CVaRUp_Neutral[i])/(CVaR_Prop[i] - CVaR_Neutral[i])
end
for i in 1:21
    CVaRUp_Compared_Neutral[i] = (CVaRUp_Prop[i] - CVaRUp_Neutral[i])/CVaRUp_Neutral[i]
end
plot_eta = plot!(X, Eta_Neutral, 
        label = "Risk-Neutral",
        color = :black,
        linewidth = s_linewidth,
        line=(:dot, 3)
)

Plots.savefig(plot_eta, results_path*"Eta_preco.png");