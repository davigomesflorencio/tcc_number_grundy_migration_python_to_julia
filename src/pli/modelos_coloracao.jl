# Modelos auxiliares de coloração própria (notebook `Original/Alg_Grundy_Number.ipynb`).
#
# No notebook original estes modelos foram usados como ponto de partida para a
# formulação de Grundy (com CP-SAT e SCIP). Aqui eles ficam como a formulação
# clássica de atribuição para o número cromático χ(G):
#
#   min  Σ_c z_c
#   s.a. Σ_c x_v^c = 1           ∀ v ∈ V
#        x_v^c ≤ z_c             ∀ v ∈ V, ∀ c ∈ C
#        x_u^c + x_v^c ≤ 1       ∀ {u,v} ∈ E, ∀ c ∈ C
#        x_v^c, z_c ∈ {0,1}

"""
    pli_numero_cromatico(g; ub=nv(g), tempo_limite=900.0, otimizador=HiGHS.Optimizer) -> Resultado

Número cromático χ(G) pela formulação de atribuição. O campo `grundy` do
`Resultado` contém χ(G). Útil como limite inferior: χ(G) ≤ Γ(G).
"""
function pli_numero_cromatico(g::SimpleGraph; ub::Integer=nv(g), tempo_limite::Real=900.0,
                              otimizador=HiGHS.Optimizer, silencioso::Bool=true)
    inicio = time()
    V = vertices(g)
    C = 1:ub
    modelo = Model(otimizador)
    silencioso && set_silent(modelo)
    set_time_limit_sec(modelo, Float64(tempo_limite))

    @variable(modelo, x[V, C], Bin)
    @variable(modelo, z[C], Bin)
    @constraint(modelo, [v in V], sum(x[v, c] for c in C) == 1)
    @constraint(modelo, [v in V, c in C], x[v, c] <= z[c])
    @constraint(modelo, [e in edges(g), c in C], x[src(e), c] + x[dst(e), c] <= 1)
    # quebra de simetria: a cor c só é usada se a cor c-1 também for
    @constraint(modelo, [c in 2:ub], z[c] <= z[c-1])
    @objective(modelo, Min, sum(z[c] for c in C))

    optimize!(modelo)
    has_values(modelo) || return Resultado(0, false, time() - inicio, Int[])
    coloracao = [something(findfirst(c -> value(x[v, c]) > 0.5, C), 0) for v in V]
    return Resultado(round(Int, objective_value(modelo)),
                     termination_status(modelo) == MOI.OPTIMAL, time() - inicio, coloracao)
end
