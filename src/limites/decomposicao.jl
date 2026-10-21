# Decomposição de vértices por grau (base do Stair Factor, da ordem de ramificação
# do BBG e da ordem externa do Bron–Kerbosch).
#
# Tradução de `ordemCrescenteGrau` + `vertex_decomposition` do TCC. O procedimento
# usa a estrutura de baldes de Batagelj–Zaversnik, mas percorre `vert` do fim para o
# início. Como a troca de posições do original só ocorre quando `u == w` (sem efeito),
# `vert` não muda e a ordem resultante é:
#   vértices em ordem decrescente de grau (empates: maior rótulo primeiro),
# e `graus_ordem[v]` = nº de vizinhos de `v` que aparecem DEPOIS de `v` na ordem.
# A tradução abaixo mantém o algoritmo original para reproduzir os resultados do TCC.

"""
    grau_maximo(g) -> Int

Grau máximo Δ(G). (Python: `calcMaxGrau`.)
"""
grau_maximo(g::SimpleGraph) = nv(g) == 0 ? 0 : Δ(g)

"""
    ordem_crescente_grau(g) -> (bin, vert, pos, graus)

Ordenação por baldes (bucket sort) dos vértices por grau. Posições são 1-based.
`bin` é indexado por grau+1. (Python: `ordemCrescenteGrau`.)
"""
function ordem_crescente_grau(g::SimpleGraph)
    n = nv(g)
    graus = degree(g)
    grau_max = n == 0 ? 0 : maximum(graus)

    bin = zeros(Int, grau_max + 1)
    for v in vertices(g)
        bin[graus[v]+1] += 1
    end
    inicio = 1
    for d in 0:grau_max
        num = bin[d+1]
        bin[d+1] = inicio
        inicio += num
    end

    pos = zeros(Int, n)
    vert = zeros(Int, n)
    for v in vertices(g)
        pos[v] = bin[graus[v]+1]
        vert[pos[v]] = v
        bin[graus[v]+1] += 1
    end

    for d in grau_max:-1:1
        bin[d+1] = bin[d]
    end
    bin[1] = 1

    return bin, vert, pos, graus
end

"""
    decomposicao_vertices(g) -> (ordem, graus_ordem)

Ordem de decomposição dos vértices e o grau residual de cada vértice no momento em
que é retirado. (Python: `vertex_decomposition`.)
"""
function decomposicao_vertices(g::SimpleGraph)
    n = nv(g)
    bin, vert, pos, graus = ordem_crescente_grau(g)
    ordem = Int[]
    graus_ordem = zeros(Int, n)

    for i in n:-1:1
        v = vert[i]
        for u in neighbors(g, v)
            du = graus[u]
            pu = pos[u]
            pw = bin[du+1]
            w = vert[pw]
            if u == w
                pos[u] = pw
                vert[pu] = w
                pos[w] = pu
                vert[pw] = u
            end
            bin[du+1] += 1
            graus[u] -= 1
        end
        graus_ordem[v] = graus[v]
        push!(ordem, v)
    end

    return ordem, graus_ordem
end
