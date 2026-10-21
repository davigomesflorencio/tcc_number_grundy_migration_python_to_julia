# Formulação de Programação Linear Inteira para o número de Grundy.
#
# Referência: http://repositorio.ufc.br/bitstream/riufc/50955/3/2020_dis_enhdrodrigues.pdf
#
#   max  Σ_c z_c
#   s.a. Σ_c x_v^c = 1                                    ∀ v ∈ V
#        z_c ≤ Σ_v x_v^c                                  ∀ c ∈ C
#        x_u^c + x_v^c ≤ z_c                              ∀ {u,v} ∈ E, ∀ c ∈ C
#        x_v^c ≥ 1 − Σ_{d<c} x_v^d − Σ_{u∈N(v)} x_u^c      ∀ v ∈ V, ∀ c ∈ C
#        x_v^c, z_c ∈ {0,1}
#
# onde C = {1, ..., UB}. As formulações do TCC diferem apenas no limite UB:
#   F_N     UB = n                   (dados_random/PLI_N.csv)
#   F_DELTA UB = Δ(G) + 1            (PLI_DELTA.csv)
#   F_STF   UB = Stair Factor        (PLI_STF.csv)
#   F_PROD  UB = 6                   (dados/PLI_6.csv, classe Stacked Book)
# (Python: `model_grundy`, resolvido com SCIP via OR-Tools.)

"""
    FORMULACOES

Formulações disponíveis: `:n`, `:delta`, `:stf`, `:prod`.
"""
const FORMULACOES = (:n, :delta, :stf, :prod)

"""
    limite_formulacao(g, formulacao) -> Int

Limite superior de cores (|C|) usado por cada formulação.
"""
function limite_formulacao(g::SimpleGraph, formulacao::Symbol)
    formulacao === :n     && return nv(g)
    formulacao === :delta && return grau_maximo(g) + 1
    formulacao === :stf   && return stair_factor(g)
    formulacao === :prod  && return 6
    throw(ArgumentError("formulação desconhecida: $formulacao (use $(FORMULACOES))"))
end

"""
    modelo_grundy(g, ub; otimizador=HiGHS.Optimizer) -> (modelo, x, z)

Constrói (sem resolver) o modelo PLI com `ub` cores.
"""
function modelo_grundy(g::SimpleGraph, ub::Integer; otimizador=HiGHS.Optimizer)
    V = vertices(g)
    C = 1:ub
    modelo = Model(otimizador)

    @variable(modelo, x[V, C], Bin)
    @variable(modelo, z[C], Bin)

    @constraint(modelo, [v in V], sum(x[v, c] for c in C) == 1)
    @constraint(modelo, [c in C], z[c] <= sum(x[v, c] for v in V))
    @constraint(modelo, [e in edges(g), c in C], x[src(e), c] + x[dst(e), c] <= z[c])
    @constraint(modelo, [v in V, c in C],
        x[v, c] >= 1 - sum(x[v, d] for d in 1:c-1; init=0) - sum(x[u, c] for u in neighbors(g, v); init=0))

    @objective(modelo, Max, sum(z[c] for c in C))
    return modelo, x, z
end

"""
    pli_grundy(g, ub; tempo_limite=900.0, otimizador=HiGHS.Optimizer, silencioso=true) -> Resultado
    pli_grundy(g, formulacao::Symbol; ...)

Resolve a formulação PLI. `ub` pode ser um inteiro ou uma das `FORMULACOES`.
Qualquer otimizador MILP do JuMP pode ser usado, por exemplo `SCIP.Optimizer`
(solver do TCC original) ou `Gurobi.Optimizer`.
"""
function pli_grundy(g::SimpleGraph, ub::Integer; tempo_limite::Real=900.0,
                    otimizador=HiGHS.Optimizer, silencioso::Bool=true)
    inicio = time()
    modelo, x, _ = modelo_grundy(g, ub; otimizador)
    silencioso && set_silent(modelo)
    set_time_limit_sec(modelo, Float64(tempo_limite))
    optimize!(modelo)

    otimo = termination_status(modelo) == MOI.OPTIMAL
    if !has_values(modelo)
        return Resultado(0, false, time() - inicio, Int[])
    end
    valor = round(Int, objective_value(modelo))
    coloracao = [something(findfirst(c -> value(x[v, c]) > 0.5, 1:ub), 0) for v in vertices(g)]
    return Resultado(valor, otimo, time() - inicio, coloracao)
end

pli_grundy(g::SimpleGraph, formulacao::Symbol; kwargs...) =
    pli_grundy(g, limite_formulacao(g, formulacao); kwargs...)
