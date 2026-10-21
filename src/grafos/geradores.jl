# Geradores de instâncias.

"""
    grafo_aleatorio(n, p; rng=Random.default_rng()) -> SimpleGraph

Grafo aleatório G(n, p): cada par {i, j} vira aresta com probabilidade `p`.
(Python: `generateGraph`.)
"""
function grafo_aleatorio(n::Integer, p::Real; rng::AbstractRNG=Random.default_rng())
    g = SimpleGraph(n)
    for i in 1:n, j in i+1:n
        rand(rng) <= p && add_edge!(g, i, j)
    end
    return g
end

"""
    grafo_stacked_book(n, m) -> SimpleGraph

Grafo da classe Stacked Book: produto cartesiano da estrela K_{1,m} com o caminho P_n.
Os vértices são numerados na ordem lexicográfica dos pares (vértice da estrela,
vértice do caminho), assim como `sorted(P1)` no Python (`generateStackedBookGraph`).
O centro da estrela é o vértice 0 da estrela.
"""
function grafo_stacked_book(n::Integer, m::Integer)
    estrela = star_graph(m + 1)   # centro = 1, folhas = 2..m+1 (nx.star_graph(m))
    caminho = path_graph(n)
    # par (h, p) -> índice lexicográfico
    indice(h, p) = (h - 1) * n + p
    g = SimpleGraph((m + 1) * n)
    for h in 1:m+1, p in 1:n
        for q in neighbors(caminho, p)
            q > p && add_edge!(g, indice(h, p), indice(h, q))
        end
        for k in neighbors(estrela, h)
            k > h && add_edge!(g, indice(h, p), indice(k, p))
        end
    end
    return g
end
