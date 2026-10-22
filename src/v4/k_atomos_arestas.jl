# BBG v4 — busca de testemunhas (k-átomos) da v3 com o limite por aresta.
#
# Mesmo algoritmo da v3 (`src/v3/k_atomos.jl`), com duas diferenças:
#
# 1. UB global = max_v b(v) de `limites_por_aresta`, nunca maior que o da v3.
# 2. Durante a busca, um vizinho u só é candidato a testemunha de cor j para o
#    vértice w se m(u→w) ≥ j (na v3 o teste era b(u) ≥ j). Isso corta candidatos
#    que só alcançariam a cor j usando o próprio w como testemunha.
#
# As demais melhorias propostas para a v4 (tabela de transposição, poda de Hall,
# simetria nas raízes, LB por busca local) ainda não estão implementadas.

"""
    bbg_v4(g; limite_superior=nothing, tempo_limite=900.0) -> Resultado

Número de Grundy exato pela busca de colorações parciais de Grundy (k-átomos),
usando os limites por aresta (`limites_por_aresta`) para o UB e para filtrar os
candidatos a testemunha.
"""
function bbg_v4(g::SimpleGraph; limite_superior::Union{Nothing,Integer}=nothing,
                tempo_limite::Real=900.0)
    inicio = time()
    nv(g) == 0 && return Resultado(0, true, 0.0, Int[])
    b, entrada = limites_por_aresta(g)
    b = min.(b, limites_por_vertice(g))
    return _grundy_por_testemunhas(g, b, entrada, limite_superior, tempo_limite, inicio)
end
