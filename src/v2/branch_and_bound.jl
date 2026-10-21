# BBG v2 — Branch and Bound melhorado para o número de Grundy.
#
# Mesma árvore de busca do BBG original (`src/exatos/branch_and_bound.jl`), com:
#
# 1. Poda válida. O original podava quando a maior cor usada ATÉ AGORA era menor que
#    o LB, mas os vértices ainda não coloridos podem receber cores maiores, então o
#    ótimo podia ser descartado. Aqui o ramo é podado quando
#        max(cor_maxima, max_{u não colorido} min(d(u)+1, UB)) ≤ LB.
# 2. Parada antecipada: se LB == UB, nenhum ramo pode melhorar e a busca termina.
# 3. LB inicial heurístico: melhor coloração gulosa entre LB-DEG e a ordem de
#    decomposição (o original começava com LB = 0).
# 4. Verificação incremental da propriedade de Grundy. Um vértice w de cor c precisa
#    de vizinhos com todas as cores 1..c-1; cada vizinho ainda sem cor fornece no
#    máximo uma cor faltante. O ramo é podado assim que
#        faltam(w) > livres(w)
#    para algum w colorido. Nas folhas livres(w) = 0, logo toda folha alcançada já é
#    uma coloração de Grundy e a verificação O(n + m) do original não é necessária.

mutable struct EstadoBBGv2
    const g::SimpleGraph{Int}
    const ordem::Vector{Int}
    const ub::Int
    const prazo::Float64
    const limite_restante::Vector{Int}   # limite_restante[i] = max_{j≥i} min(d(v_j)+1, UB)
    const coloracao::Vector{Int}
    const contagem::Matrix{Int}          # contagem[k, w] = nº de vizinhos de w com cor k
    const faltam::Vector{Int}            # cores 1..c(w)-1 ausentes na vizinhança de w
    const livres::Vector{Int}            # vizinhos de w ainda sem cor
    lb::Int
    melhor_coloracao::Vector{Int}
    nos::Int
    parar::Bool
end

"""
    bbg_v2(g, limite_superior; tempo_limite=900.0, lb_heuristico=true) -> Resultado

Branch and Bound melhorado (ver o cabeçalho de `src/v2/branch_and_bound.jl`).
Recebe os mesmos limites do `bbg` original (Δ+1, Stair Factor, 6) e, ao contrário
dele, é exato: devolve o maior número de cores de uma coloração de Grundy que usa
no máximo `limite_superior` cores (ou o LB heurístico, se este já for maior).
Com UB = Δ+1 ou Stair Factor, que são limites superiores válidos, o resultado é o
número de Grundy.
"""
function bbg_v2(g::SimpleGraph, limite_superior::Integer; tempo_limite::Real=900.0,
                lb_heuristico::Bool=true)
    inicio = time()
    n = nv(g)
    ub = Int(limite_superior)
    ordem, _ = decomposicao_vertices(g)

    limite_restante = zeros(Int, n + 1)
    for i in n:-1:1
        limite_restante[i] = max(limite_restante[i+1], min(degree(g, ordem[i]) + 1, ub))
    end

    estado = EstadoBBGv2(g, ordem, ub, inicio + tempo_limite, limite_restante,
                         fill(SEM_COR, n), zeros(Int, max(ub, 1), n), zeros(Int, n),
                         degree(g), 0, Int[], 0, false)

    if lb_heuristico
        for (k, c) in (lb_deg(g), coloracao_gulosa(g, ordem))
            k > estado.lb && ((estado.lb, estado.melhor_coloracao) = (k, c))
        end
    end

    otimo = true
    try
        estado.lb < ub && _bbg_v2!(estado, 1, 0)
    catch e
        e isa TempoEsgotado || rethrow()
        otimo = false
    end
    return Resultado(estado.lb, otimo, time() - inicio, estado.melhor_coloracao)
end

function _bbg_v2!(estado::EstadoBBGv2, profundidade::Int, cor_maxima::Int)
    estado.nos += 1
    (estado.nos & 0x3ff) == 0 && _verificar_prazo(estado.prazo)

    # (1) poda por limite válido
    max(cor_maxima, estado.limite_restante[profundidade]) <= estado.lb && return

    if profundidade > length(estado.ordem)
        # (4) garante que a folha é uma coloração de Grundy
        estado.lb = cor_maxima
        estado.melhor_coloracao = copy(estado.coloracao)
        estado.parar = estado.lb >= estado.ub          # (2) parada antecipada
        return
    end

    g, coloracao, contagem = estado.g, estado.coloracao, estado.contagem
    v = estado.ordem[profundidade]
    limite = min(degree(g, v) + 1, estado.ub)

    for cor in limite:-1:1
        contagem[cor, v] > 0 && continue                # cor usada por um vizinho
        if _atribuir!(estado, v, cor)
            _bbg_v2!(estado, profundidade + 1, max(cor_maxima, cor))
        end
        _desfazer!(estado, v, cor)
        estado.parar && return
        # (1) com a cor maior já testada, as menores só podem alcançar o mesmo limite
        max(cor_maxima, estado.limite_restante[profundidade+1], cor - 1) <= estado.lb && return
    end
end

# Colore `v` com `cor`, atualiza os contadores e devolve `false` se algum vértice
# afetado (v ou um vizinho colorido) não puder mais completar suas cores menores.
function _atribuir!(estado::EstadoBBGv2, v::Int, cor::Int)
    g, coloracao, contagem, faltam, livres = estado.g, estado.coloracao, estado.contagem,
                                              estado.faltam, estado.livres
    coloracao[v] = cor
    viavel = true
    for u in neighbors(g, v)
        livres[u] -= 1
        contagem[cor, u] += 1
        cu = coloracao[u]
        if cu != SEM_COR
            cor < cu && contagem[cor, u] == 1 && (faltam[u] -= 1)
            faltam[u] > livres[u] && (viavel = false)
        end
    end
    f = 0
    for k in 1:cor-1
        contagem[k, v] == 0 && (f += 1)
    end
    faltam[v] = f
    return viavel && f <= livres[v]
end

function _desfazer!(estado::EstadoBBGv2, v::Int, cor::Int)
    g, coloracao, contagem, faltam, livres = estado.g, estado.coloracao, estado.contagem,
                                              estado.faltam, estado.livres
    for u in neighbors(g, v)
        livres[u] += 1
        contagem[cor, u] -= 1
        cu = coloracao[u]
        cu != SEM_COR && cor < cu && contagem[cor, u] == 0 && (faltam[u] += 1)
    end
    faltam[v] = 0
    coloracao[v] = SEM_COR
end
