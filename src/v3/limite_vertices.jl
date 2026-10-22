# Limite superior por vértice para a cor que cada vértice pode receber em uma
# coloração de Grundy.
#
# Se v recebe a cor c, seus vizinhos precisam ter as cores 1..c-1, e um vizinho de
# cor j precisa, por sua vez, poder receber a cor j. Partindo de b(v) = d(v) + 1,
# itera-se até o ponto fixo:
#
#   b(v) = 1 + max{ t ≤ d(v) : existem vizinhos distintos u_1..u_t com b(u_i) ≥ i }
#
# Com os limites dos vizinhos em ordem decrescente β_1 ≥ β_2 ≥ ..., o valor t é
# viável se e só se β_i ≥ t - i + 1 para i = 1..t (as cores maiores vão para os
# vizinhos de maior limite). Por indução sobre as iterações, toda coloração de
# Grundy satisfaz cor(v) ≤ b(v); logo Γ(G) ≤ max_v b(v) — diferente do Stair Factor
# do TCC, este é sempre um limite válido.

"""
    limites_por_vertice(g) -> Vector{Int}

`b[v]` = limite superior para a cor de `v` em qualquer coloração de Grundy.
"""
function limites_por_vertice(g::SimpleGraph)
    b = degree(g) .+ 1
    vizinhos_b = Int[]
    mudou = true
    while mudou
        mudou = false
        for v in vertices(g)
            empty!(vizinhos_b)
            for u in neighbors(g, v)
                push!(vizinhos_b, b[u])
            end
            sort!(vizinhos_b; rev=true)
            novo = 1 + _maior_t_viavel(vizinhos_b)
            if novo < b[v]
                b[v] = novo
                mudou = true
            end
        end
    end
    return b
end

# maior t com β[i] ≥ t - i + 1 para todo i ≤ t  ⇔  t ≤ min_{i≤t} (β[i] + i - 1)
function _maior_t_viavel(β::Vector{Int})
    melhor = 0
    minimo = typemax(Int)
    for t in eachindex(β)
        minimo = min(minimo, β[t] + t - 1)
        t <= minimo && (melhor = t)
    end
    return melhor
end

"""
    limite_superior_vertices(g) -> Int

Limite superior válido para Γ(G): `max_v b(v)`.
"""
limite_superior_vertices(g::SimpleGraph) = nv(g) == 0 ? 0 : maximum(limites_por_vertice(g))
