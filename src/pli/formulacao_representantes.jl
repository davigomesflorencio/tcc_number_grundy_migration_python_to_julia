# Formulação de representantes para o número de Grundy (padrão, não conexo).
#
# Referência: Carvalho, Melo, Santos, Toso e Resende — formulação de representantes
# para a coloração de Grundy (SBPO 2023; versão conexa em arXiv:2411.14533, Seção 2.2).
# A versão abaixo é a versão padrão: sem o índice de tempo t e sem as restrições de
# conectividade.
#
# Notação: N̄(u) = V \ N[u] (antivizinhança), N̄[u] = N̄(u) ∪ {u}. A cor de cada classe
# é representada pelo vértice de menor índice da classe (u ≤ v).
#
# Variáveis:
#   x[u,v] = 1 se u representa v      (u ≤ v, v ∈ N̄[u]); x[u,u] = 1 ⇔ u é representante
#   y[p,u] = 1 se p e u são representantes e a classe de p vem antes da classe de u
#
#   max  Σ_u x[u,u]
#   s.a. Σ_{u ∈ N̄[v], u ≤ v} x[u,v] = 1                         ∀ v ∈ V
#        x[u,v] + x[u,w] ≤ x[u,u]                               ∀ u, {v,w} ∈ E, v,w ∈ N̄(u), u < v,w
#        x[u,v] ≤ x[u,u]                                        ∀ u < v, v ∈ N̄(u)
#        x[u,v] ≤ Σ_{w ∈ N(v) ∩ N̄[p], p ≤ w} x[p,w] + 1 − y[p,u]  ∀ p ≠ u, v ∈ N̄[u], u ≤ v   (Grundy)
#        y[p,u] + y[u,p] ≥ x[p,p] + x[u,u] − 1                  ∀ p < u
#        y[p,u] + y[u,p] ≤ x[p,p],  y[p,u] ≤ x[u,u]             ∀ p ≠ u
#        y[p,u] + y[u,q] − y[p,q] ≤ 1                           ∀ p, u, q distintos   (transitividade)
#        Σ_u x[u,u] ≤ UB
#        x, y ∈ {0,1}
#
# A transitividade não aparece no artigo, mas sem ela um ciclo A → B → C → A entre as
# classes seria aceito mesmo sem existir ordem de Grundy válida; com ela, y é uma
# ordem total sobre os representantes.

"""
    modelo_representantes(g, ub; otimizador=HiGHS.Optimizer) -> (modelo, x, y)

Constrói (sem resolver) a formulação de representantes com no máximo `ub` cores.
`x` e `y` são dicionários indexados por pares `(u, v)`.
"""
function modelo_representantes(g::SimpleGraph, ub::Integer; otimizador=HiGHS.Optimizer)
    V = vertices(g)
    adj = [Set(neighbors(g, v)) for v in V]
    anti(u, v) = u == v || !(v in adj[u])          # v ∈ N̄[u]
    pares = [(u, v) for v in V for u in V if u <= v && anti(u, v)]
    modelo = Model(otimizador)

    @variable(modelo, x[pares], Bin)
    @variable(modelo, y[[(p, u) for p in V for u in V if p != u]], Bin)

    # cada vértice tem exatamente um representante
    @constraint(modelo, [v in V], sum(x[(u, v)] for u in 1:v if anti(u, v)) == 1)
    # vértices adjacentes não compartilham representante
    for u in V, e in edges(g)
        v, w = src(e), dst(e)
        if u < v && u < w && anti(u, v) && anti(u, w)
            @constraint(modelo, x[(u, v)] + x[(u, w)] <= x[(u, u)])
        end
    end
    # só representantes representam
    @constraint(modelo, [(u, v) in pares; u < v], x[(u, v)] <= x[(u, u)])
    # propriedade de Grundy: se p vem antes de u, todo vértice da classe de u
    # tem um vizinho na classe de p
    for p in V, (u, v) in pares
        p == u && continue
        @constraint(modelo, x[(u, v)] <=
            sum(x[(p, w)] for w in adj[v] if p <= w && anti(p, w); init=0) + 1 - y[(p, u)])
    end
    # ordem total entre os representantes
    for p in V, u in V
        p == u && continue
        p < u && @constraint(modelo, y[(p, u)] + y[(u, p)] >= x[(p, p)] + x[(u, u)] - 1)
        @constraint(modelo, y[(p, u)] + y[(u, p)] <= x[(p, p)])
        @constraint(modelo, y[(p, u)] <= x[(u, u)])
        for q in V
            (q == p || q == u) && continue
            @constraint(modelo, y[(p, u)] + y[(u, q)] - y[(p, q)] <= 1)
        end
    end
    @constraint(modelo, sum(x[(u, u)] for u in V) <= ub)

    @objective(modelo, Max, sum(x[(u, u)] for u in V))
    return modelo, x, y
end

"""
    pli_representantes(g, ub; tempo_limite=900.0, otimizador=HiGHS.Optimizer, silencioso=true) -> Resultado
    pli_representantes(g, formulacao::Symbol; ...)
    pli_representantes(g; ...)                       # ub = Stair Factor

Resolve a formulação de representantes. `ub` pode ser um inteiro ou uma das
`FORMULACOES` (`:n`, `:delta`, `:stf`, `:prod`).
"""
function pli_representantes(g::SimpleGraph, ub::Integer; tempo_limite::Real=900.0,
                            otimizador=HiGHS.Optimizer, silencioso::Bool=true)
    inicio = time()
    modelo, x, y = modelo_representantes(g, ub; otimizador)
    silencioso && set_silent(modelo)
    set_time_limit_sec(modelo, Float64(tempo_limite))
    optimize!(modelo)

    otimo = termination_status(modelo) == MOI.OPTIMAL
    if !has_values(modelo)
        return Resultado(0, false, time() - inicio, Int[])
    end
    valor = round(Int, objective_value(modelo))

    # cor da classe de u = 1 + número de representantes que vêm antes de u
    V = vertices(g)
    representa(u, v) = (u == v || !has_edge(g, u, v)) && value(x[(u, v)]) > 0.5
    rep(v) = something(findfirst(u -> representa(u, v), 1:v), 0)
    reps = [u for u in V if value(x[(u, u)]) > 0.5]
    cor = Dict(u => 1 + count(p -> p != u && value(y[(p, u)]) > 0.5, reps) for u in reps)
    coloracao = [get(cor, rep(v), 0) for v in V]
    return Resultado(valor, otimo, time() - inicio, coloracao)
end

pli_representantes(g::SimpleGraph, formulacao::Symbol; kwargs...) =
    pli_representantes(g, limite_formulacao(g, formulacao); kwargs...)

pli_representantes(g::SimpleGraph; kwargs...) = pli_representantes(g, :stf; kwargs...)
