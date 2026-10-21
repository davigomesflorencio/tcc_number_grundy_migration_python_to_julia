# BBG — Algoritmo Branch and Bound para o número de Grundy.
#
# Os vértices são coloridos na ordem da decomposição de vértices. Para cada vértice
# `v` testam-se as cores 1..min(grau(v)+1, UB) não usadas pelos vizinhos já
# coloridos, da maior para a menor. Um ramo é podado quando a maior cor usada até o
# momento é menor que o melhor valor conhecido (LB). Nas folhas verifica-se se a
# coloração é de Grundy. (Python: `BAB` + `BABNumberGrundy`.)

mutable struct EstadoBBG
    const g::SimpleGraph{Int}
    const ordem::Vector{Int}
    const ub::Int
    const prazo::Float64
    lb::Int
    melhor_coloracao::Vector{Int}
    nos::Int
end

"""
    bbg(g, limite_superior; tempo_limite=900.0) -> Resultado

Número de Grundy via Branch and Bound. `limite_superior` limita as cores testadas
(no TCC: Δ+1 → BB_DELTA, Stair Factor → BB_STF, 6 → BB_6).
"""
function bbg(g::SimpleGraph, limite_superior::Integer; tempo_limite::Real=900.0)
    inicio = time()
    ordem, _ = decomposicao_vertices(g)
    estado = EstadoBBG(g, ordem, limite_superior, inicio + tempo_limite, 0, Int[], 0)
    coloracao = fill(SEM_COR, nv(g))
    otimo = true
    try
        _bbg!(estado, coloracao, 1, 0)
    catch e
        e isa TempoEsgotado || rethrow()
        otimo = false
    end
    return Resultado(estado.lb, otimo, time() - inicio, estado.melhor_coloracao)
end

function _bbg!(estado::EstadoBBG, coloracao::Vector{Int}, profundidade::Int, cor_maxima::Int)
    g = estado.g
    estado.nos += 1
    (estado.nos & 0x3ff) == 0 && _verificar_prazo(estado.prazo)

    if profundidade > length(estado.ordem)
        if cor_maxima > estado.lb && eh_coloracao_grundy(g, coloracao)
            estado.lb = cor_maxima
            estado.melhor_coloracao = copy(coloracao)
        end
        return
    end

    v = estado.ordem[profundidade]
    limite = min(degree(g, v) + 1, estado.ub)
    proibidas = falses(limite)
    for u in neighbors(g, v)
        c = coloracao[u]
        1 <= c <= limite && (proibidas[c] = true)
    end

    for cor in limite:-1:1
        proibidas[cor] && continue
        nova_max = max(cor_maxima, cor)
        if nova_max >= estado.lb
            coloracao[v] = cor
            _bbg!(estado, coloracao, profundidade + 1, nova_max)
            coloracao[v] = SEM_COR
        end
    end
end
