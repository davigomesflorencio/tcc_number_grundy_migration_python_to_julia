# Limite superior Stair Factor.

"""
    stair_factor(g) -> Int

Limite superior para o número de Grundy:
`min_i ( i + graus_ordem[v_i] )` sobre a ordem de decomposição `v_1, ..., v_n`.
(Python: `stairfactor`.)
"""
function stair_factor(g::SimpleGraph)
    ordem, graus_ordem = decomposicao_vertices(g)
    isempty(ordem) && return 0
    return minimum(i + graus_ordem[v] for (i, v) in enumerate(ordem))
end
