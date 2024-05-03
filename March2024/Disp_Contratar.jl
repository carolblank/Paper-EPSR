
# ----------------------------- Upside -------------------------------

range_P = collect(0:10:200)

x_otimo_disp_upside = zeros(length(range_P))
CVaR_Prop = zeros(length(range_P))
CVaRUp_Prop = zeros(length(range_P))
CVaR_Prop_Test = zeros(length(range_P))
CVaRUp_Prop_Test = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel

        λ = 1.0;

        global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
        Rdisp_otimo1, Rquant_otimo1, R_otimo1, R_total_otimo1, x_otimo1, y_otimo1 = port_alloc(params);

        cvar    = Statistics.mean(sort(R_total_otimo1)[1:Int(round((1 - α) * S)), :]);
        avg_adj = mean(R_total_otimo1);
        L = cvar/pu_Money;
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


        S_Test         = 800
        R_teste_upside = zeros(S_Test)

        for s in 1:S_Test
                R_teste_upside[s] = sum(sum(gu[s+S, t, i] * GF[i] * PLD[s+S, t] * h[t] for i in 1:nI) 
                                + sum((P[j, t] - PLD[s+S, t]) * Q[j, t] * x_novo[j] * h[t] for j in 1:J) for t in 1:T)./1e6;
        end

        CVaR_Prop_Test[iter]       = Statistics.mean(sort(R_teste_upside)[1:Int(round((1 - α) * S_Test)), :]);
        CVaRUp_Prop_Test[iter]       = Statistics.mean(sort(R_teste_upside)[Int(round((α) * S_Test))+1:S_Test,:]);

end

Y = x_otimo_disp_upside
X = range_P

p = plot(X, Y, 
        label = "Proposed approach",
        xlabel = "Contract Price (A+1)",
        ylabel = "Contracting Level",
        color=:blue,
        linewidth = s_linewidth
)




# ----------------------------- Risk-Reward -------------------------------


x_otimo_disp_RR = zeros(length(range_P))
CVaR_RR = zeros(length(range_P))
CVaRUp_RR = zeros(length(range_P))
CVaR_RR_Test = zeros(length(range_P))
CVaRUp_RR_Test = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel
        S_Test         = 800
        R_teste_RR = zeros(S_Test)

        λ_ = 0.5

        global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
        Rdisp_RR, Rquant_RR, RPort_RR, RPortTot_RR, xOptimal_RR, yOptimal_RR = port_alloc(params);
    
        x_otimo_disp_RR[iter] = xOptimal_RR[1]
        CVaR_RR[iter]             = Statistics.mean(sort(RPortTot_RR)[1:Int(round((1 - α) * S)), :])/pu_Money;
        CVaRUp_RR[iter]           = Statistics.mean(sort(RPortTot_RR)[Int(round((α) * S))+1:S, :])/pu_Money;

        for s in 1:S_Test
                R_teste_RR[s] = sum(sum(gu[s+S, t, i] * GF[i] * PLD[s+S, t] * h[t] for i in 1:nI) 
                                + sum((P[j, t] - PLD[s+S, t]) * Q[j, t] * xOptimal_RR[j] * h[t] for j in 1:J) for t in 1:T);
        end
        
        CVaR_RR_Test[iter]             = Statistics.mean(sort(R_teste_RR)[1:Int(round((1 - α) * S_Test)), :])/pu_Money;
        CVaRUp_RR_Test[iter]           = Statistics.mean(sort(R_teste_RR)[Int(round((α) * S_Test))+1:S_Test, :])/pu_Money;

end

X = range_P
Y = x_otimo_disp_RR

p = plot!(X, Y, 
        label = "Risk-Averse (λ = 0.5)",
        color = :gray,
        linewidth = s_linewidth,
        ls=:dash
)


# ----------------------------- Risk-neutral -------------------------------

x_otimo_disp_Neutral = zeros(length(range_P))
CVaR_Neutral = zeros(length(range_P))
CVaRUp_Neutral = zeros(length(range_P))
CVaR_Neutral_Test = zeros(length(range_P))
CVaRUp_Neutral_Test = zeros(length(range_P))

for (iter,p_variavel) in enumerate(range_P)

        P[1,:]     .= p_variavel

        λ_ = 0.0
        S_Test         = 800
        R_teste_Neutral = zeros(S_Test)

        global params = ParametrosContratacao(gu, GF, h, Q, PLD, α, T, λ_, flag_restricao, Rmin, P, S, q, nI, J, C, tx_m, xmin, xmax, ymin, ymax, Solv);
        Rdisp_Neutral, Rquant_Neutral, RPort_Neutral, RPortTot_Neutral, xOptimal_Neutral, yOptimal_Neutral = port_alloc(params);
    
        x_otimo_disp_Neutral[iter] = xOptimal_Neutral[1]
        CVaR_Neutral[iter]             = Statistics.mean(sort(RPortTot_Neutral)[1:Int(round((1 - α) * S)), :])/pu_Money;
        CVaRUp_Neutral[iter]           = Statistics.mean(sort(RPortTot_Neutral)[Int(round((α) * S))+1:S, :])/pu_Money;

        for s in 1:S_Test
                R_teste_Neutral[s] = sum(sum(gu[s+S, t, i] * GF[i] * PLD[s+S, t] * h[t] for i in 1:nI) 
                                + sum((P[j, t] - PLD[s+S, t]) * Q[j, t] * xOptimal_Neutral[j] * h[t] for j in 1:J) for t in 1:T);
        end
        
        CVaR_Neutral_Test[iter]             = Statistics.mean(sort(R_teste_Neutral)[1:Int(round((1 - α) * S_Test)), :])/pu_Money;
        CVaRUp_Neutral_Test[iter]           = Statistics.mean(sort(R_teste_Neutral)[Int(round((α) * S_Test))+1:S_Test, :])/pu_Money;

end


X = range_P
Y = x_otimo_disp_Neutral

p = plot!(X,Y, 
        label = "Risk-Neutral",
        color = :black,
        linewidth = s_linewidth,
        line=(:dot, 3)
)

display(p)

Plots.savefig(plot_eta, results_path*"Eta_preco.png");

Eta_RR = zeros(21)
CVaRUp_Compared_RR = zeros(21)
Eta_RR = (CVaRUp_Prop - CVaRUp_RR) ./ (CVaR_Prop - CVaR_RR)
for i in 1:21
    CVaRUp_Compared_RR[i] = (CVaRUp_Prop[i] - CVaRUp_RR[i])/CVaRUp_RR[i]
end
plot_eta = plot(X[10:end], Eta_RR[10:end], 
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
plot_eta = plot!(X[10:end], Eta_Neutral[10:end], 
        label = "Risk-Neutral",
        color = :black,
        linewidth = s_linewidth,
        line=(:dot, 3)
)





Eta_RR_Test = zeros(10)
for i in 10:19
    if CVaRUp_Prop_Test[i] - CVaRUp_RR_Test[i] < 10e-5
        Eta_RR_Test[i-9] = 0
    else
        Eta_RR_Test[i-9] = (CVaRUp_Prop_Test[i] - CVaRUp_RR_Test[i])/(CVaR_Prop_Test[i] - CVaR_RR_Test[i])
    end
end
a = CVaRUp_Prop_Test[10:end-2] - CVaRUp_RR_Test[10:end-2] 
b = CVaR_Prop_Test[10:end-2] - CVaR_RR_Test[10:end-2]

Eta_Neutral_Test = zeros(10)
for i in 10:19
    Eta_Neutral_Test[i-9] = (CVaRUp_Prop_Test[i] - CVaRUp_Neutral_Test[i])/(CVaR_Prop_Test[i] - CVaR_Neutral_Test[i])
end
c = CVaRUp_Prop_Test[10:end-2] - CVaRUp_Neutral_Test[10:end-2]
d = CVaR_Prop_Test[10:end-2] - CVaR_Neutral_Test[10:end-2]

column_titles = ["90" ; "100" ; "110" ; "120" ; "130" ; "140" ; "150" ; "160" ; "170"; "180"]
qMat    = [a b Eta_RR_Test c d Eta_Neutral_Test]';
q_df    = DataFrame(qMat, :auto);
q_df = DataFrame(qMat, Symbol.(column_titles))

row_titles = ["CVaU Prop RR"; "CVaR Prop RR"; "Eta Prop RR"; "CVaU Prop Neutral"; "CVaR Prop Neutral"; "Eta Prop Neutral"]
insertcols!(q_df, 1, "Métrica" => row_titles)


Eta_Neutral_Test = zeros(21)
for i in 1:21
    Eta_Neutral_Test[i] = (CVaRUp_Prop_Test[i] - CVaRUp_Neutral_Test[i])/(CVaR_Prop_Test[i] - CVaR_Neutral_Test[i])
end