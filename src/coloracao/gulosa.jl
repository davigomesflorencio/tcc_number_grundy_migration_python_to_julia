# Coloração gulosa e verificação de coloração de Grundy.
# Cores são inteiros 1, 2, ...; `0` indica vértice ainda não colorido
# (no Python era `-1`).

const SEM_COR = 0

"""
    cor_minima(g, v, coloracao) -> Int

Menor cor positiva não usada pelos vizinhos já coloridos de `v`. (Python: `corMinima`.)
"""
function cor_minima(g::SimpleGraph, v::Integer, coloracao::AbstractVector{<:Integer})
    usadas = falses(degree(g, v) + 1)
    for u in neighbors(g, v)
        c = coloracao[u]
        1 <= c <= length(usadas) && (usadas[c] = true)
    end
    return findfirst(!, usadas)
end

"""
    coloracao_gulosa(g, ordem) -> (max_cor, coloracao)

Algoritmo Guloso de Coloração: percorre os vértices em `ordem` atribuindo sempre a
menor cor possível. (Python: `coloracaoGulosa`.)
"""
function coloracao_gulosa(g::SimpleGraph, ordem::AbstractVector{<:Integer}=collect(vertices(g)))
    coloracao = fill(SEM_COR, nv(g))
    for v in ordem
        coloracao[v] = cor_minima(g, v, coloracao)
    end
    return (isempty(coloracao) ? 0 : maximum(coloracao)), coloracao
end

"""
    eh_coloracao_grundy(g, coloracao) -> Bool

Verifica se todo vértice de cor `c` possui vizinhos com todas as cores `1..c-1`.
(Python: `isGrundyColoring`.)
"""
function eh_coloracao_grundy(g::SimpleGraph, coloracao::AbstractVector{<:Integer})
    for v in vertices(g)
        cor = coloracao[v]
        vistas = falses(max(cor - 1, 0))
        for u in neighbors(g, v)
            c = coloracao[u]
            1 <= c < cor && (vistas[c] = true)
        end
        all(vistas) || return false
    end
    return true
end
