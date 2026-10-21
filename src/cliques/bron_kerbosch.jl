# Algoritmo de Bron–Kerbosch para enumeração de cliques maximais.
#
# Duas versões:
# - `cliques_maximais(g)`: versão com pivô e ordem externa dada pela decomposição de
#   vértices (Python: `findCliques` + `findCliques3`).
# - `conjuntos_independentes_maximais(adj, S, f)`: versão em bitmask usada pela
#   Programação Dinâmica; enumera os cliques maximais do complemento de G[S], isto é,
#   os conjuntos independentes maximais de G[S].

"""
    cliques_maximais(g) -> Vector{Vector{Int}}

Todos os cliques maximais de `g`.
"""
function cliques_maximais(g::SimpleGraph)
    cliques = Vector{Int}[]
    vizinhos = [Set(neighbors(g, v)) for v in vertices(g)]
    P = Set(vertices(g))
    X = Set{Int}()
    ordem, _ = decomposicao_vertices(g)
    for v in Iterators.reverse(ordem)
        _bron_kerbosch!(cliques, vizinhos, intersect(P, vizinhos[v]), [v], intersect(X, vizinhos[v]))
        delete!(P, v)
        push!(X, v)
    end
    return cliques
end

function _bron_kerbosch!(cliques, vizinhos, P::Set{Int}, R::Vector{Int}, X::Set{Int})
    if isempty(P) && isempty(X)
        push!(cliques, sort(R))
        return
    end
    # pivô: vértice de P ∪ X com mais vizinhos em P
    pivo = argmax(u -> length(intersect(P, vizinhos[u])), union(P, X))
    for v in collect(setdiff(P, vizinhos[pivo]))
        _bron_kerbosch!(cliques, vizinhos, intersect(P, vizinhos[v]), [R; v], intersect(X, vizinhos[v]))
        delete!(P, v)
        push!(X, v)
    end
end

"""
    mascaras_adjacencia(g) -> Vector{UInt64}

`adj[v]` tem o bit `v-1` ligado para cada vizinho de `v`. Requer `nv(g) <= 64`.
"""
function mascaras_adjacencia(g::SimpleGraph)
    nv(g) <= 64 || throw(ArgumentError("representação em bitmask requer até 64 vértices"))
    adj = zeros(UInt64, nv(g))
    for v in vertices(g), u in neighbors(g, v)
        adj[v] |= UInt64(1) << (u - 1)
    end
    return adj
end

"""
    conjuntos_independentes_maximais(f, adj, S)

Chama `f(I)` para cada conjunto independente maximal `I` (bitmask) do subgrafo
induzido por `S` (bitmask), com `adj` vindo de `mascaras_adjacencia`.
"""
function conjuntos_independentes_maximais(f::F, adj::Vector{UInt64}, S::UInt64) where {F}
    S == 0 && return nothing
    _bk_independentes(f, adj, UInt64(0), S, UInt64(0))
    return nothing
end

# Bron–Kerbosch com pivô no complemento: os "vizinhos" de v são os não-adjacentes.
@inline _nao_vizinhos(adj, v) = ~adj[v] & ~(UInt64(1) << (v - 1))

function _bk_independentes(f::F, adj, R::UInt64, P::UInt64, X::UInt64) where {F}
    if P == 0
        X == 0 && f(R)
        return
    end
    # pivô: vértice de P ∪ X que maximiza |P ∩ N̄(u)|
    PX = P | X
    pivo, melhor = 0, -1
    while PX != 0
        u = trailing_zeros(PX) + 1
        k = count_ones(P & _nao_vizinhos(adj, u))
        k > melhor && ((pivo, melhor) = (u, k))
        PX &= PX - 1
    end
    candidatos = P & ~_nao_vizinhos(adj, pivo)
    while candidatos != 0
        v = trailing_zeros(candidatos) + 1
        bit = UInt64(1) << (v - 1)
        nviz = _nao_vizinhos(adj, v)
        _bk_independentes(f, adj, R | bit, P & nviz, X & nviz)
        P &= ~bit
        X |= bit
        candidatos &= candidatos - 1
    end
end
