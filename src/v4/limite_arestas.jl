# Limite superior por aresta ("non-backtracking") para a cor dos vértices em uma
# coloração de Grundy.
#
# Para cada aresta orientada u→v, m(u→v) limita a cor de u quando v tem cor MAIOR
# que a de u — é o caso em que u é testemunha de v. Nessa situação v não pode ser
# testemunha de u, então só os demais vizinhos de u contam:
#
#   m(u→v) = 1 + T( { m(x→u) : x ∈ N(u) \ {v} } )
#   b(v)   = 1 + T( { m(u→v) : u ∈ N(v) } )
#
# onde T(β) = max{ t : com β em ordem decrescente, β_i ≥ t - i + 1 para i = 1..t }
# (mesma regra de `limites_por_vertice`). Parte-se de m(u→v) = d(u) — u tem d(u)-1
# outros vizinhos — e itera-se até o ponto fixo.
#
# Validade, por indução na cor de u: se cor(u) < cor(v), as testemunhas de u são
# vizinhos x ≠ v com cor(x) < cor(u), logo cor(x) ≤ m(x→u) e cor(u) ≤ m(u→v).
# Portanto cor(v) ≤ b(v) e Γ(G) ≤ max_v b(v).
#
# Em relação à v3 (`limites_por_vertice`), o ganho é não deixar uma testemunha
# "usar" o próprio vértice que ela testemunha. Ex.: na estrela K_{1,m} a v3 dá 3 e
# este limite dá 2 (= Γ); na Stacked Book K_{1,m} □ P_n dá 6 em vez de 7.

"""
    limites_por_aresta(g) -> (b, entrada)

- `b[v]`: limite superior para a cor de `v` em qualquer coloração de Grundy;
- `entrada[v][i]`: limite m(u→v) para a cor do i-ésimo vizinho `u` de `v` quando
  `u` é testemunha de `v` (cor de u < cor de v).

Custo por iteração: O(Σ_v d(v)² log d(v)).
"""
function limites_por_aresta(g::SimpleGraph)
    n = nv(g)
    # entrada[v][i] = m(u→v), u = neighbors(g, v)[i]
    entrada = [degree(g)[neighbors(g, v)] for v in vertices(g)]
    # posicao[v][i] = índice de v em neighbors(g, u), u = neighbors(g, v)[i]
    posicao = [[searchsortedfirst(neighbors(g, u), v) for u in neighbors(g, v)] for v in vertices(g)]

    buffer = Int[]
    mudou = true
    while mudou
        mudou = false
        for u in vertices(g)
            viz = neighbors(g, u)
            # nova mensagem m(u→v) para cada vizinho v, excluindo v da lista
            for (i, v) in enumerate(viz)
                empty!(buffer)
                for (j, x) in enumerate(viz)
                    j != i && push!(buffer, entrada[u][j])
                end
                sort!(buffer; rev=true)
                novo = 1 + _maior_t_viavel(buffer)
                p = posicao[u][i]                  # u é o p-ésimo vizinho de v
                if novo < entrada[v][p]
                    entrada[v][p] = novo
                    mudou = true
                end
            end
        end
    end

    b = Vector{Int}(undef, n)
    for v in vertices(g)
        buffer = sort(entrada[v]; rev=true)
        b[v] = 1 + _maior_t_viavel(buffer)
    end
    return b, entrada
end

"""
    limite_superior_arestas(g) -> Int

Limite superior válido para Γ(G): `max_v b(v)` com `b` de `limites_por_aresta`.
"""
limite_superior_arestas(g::SimpleGraph) = nv(g) == 0 ? 0 : maximum(limites_por_aresta(g)[1])
