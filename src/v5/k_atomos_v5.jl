# BBG v5 — busca de k-átomos da v4 com tabela de transposição e poda de Hall.
#
# Mesma busca da v4 (`src/v4/k_atomos_arestas.jl`): limite por aresta para o UB e
# para filtrar candidatos a testemunha. Acrescenta:
#
# 1. Tabela de transposição (memoização de estados que falharam). A mesma
#    coloração parcial pode ser alcançada por caminhos diferentes da busca (colorir
#    u e depois x, ou x e depois u, quando atendem demandas diferentes). Cada estado
#    tem um hash Zobrist de 128 bits, atualizado em O(1) a cada (des)coloração. Se a
#    busca a partir de um estado falha, o hash é guardado e o estado não é explorado
#    de novo. Isso é válido porque a falha depende só do estado, de k e das raízes
#    proibidas — e as proibições só aumentam durante um teste. A tabela é zerada a
#    cada novo k e tem tamanho máximo (`max_transposicoes`); ao atingir o máximo,
#    novas falhas deixam de ser guardadas (a busca continua correta). Só são
#    guardados estados cuja subárvore teve pelo menos `min_subarvore` nós: guardar
#    folhas baratas custa mais do que refazê-las.
#
# 2. Poda de Hall. Um vértice colorido w precisa receber, de vizinhos distintos
#    ainda sem cor, cada cor que lhe falta. Monta-se o grafo bipartido
#    (cores faltantes de w) × (vizinhos sem cor de w que podem receber a cor) e
#    procura-se um emparelhamento que cubra todas as cores faltantes (algoritmo de
#    caminhos aumentantes de Kuhn, O(k · arestas) por vértice). Se não existe, o ramo
#    é podado. Substitui as podas da v3 "faltam > livres" e "cor sem candidato",
#    que são casos particulares do teorema de Hall. O emparelhamento só é calculado
#    quando alguma cor faltante tem menos de `faltam` candidatos: caso contrário a
#    condição de Hall vale automaticamente.
#
# A chance de colisão de dois estados diferentes em 128 bits é desprezível
# (~ q² / 2^129 para q estados guardados).
#
# Medições (random_1_20_0.5): v4 18,1 M nós / 9,7 s; só transposição 15,4 M nós
# (−15%) / 9,2 s; só Hall 17,6 M nós (−2,5%) / 12,9 s. A poda de Hall quase não
# reduz a árvore porque a escolha da demanda mais restrita (MRV) já detecta o mesmo
# conflito um nó depois, e o emparelhamento custa mais que esse nó. Por isso a
# poda de Hall vem DESLIGADA por padrão (`hall=false`); a transposição, ligada.

mutable struct EstadoV5
    const g::SimpleGraph{Int}
    const k::Int
    const b::Vector{Int}
    const entrada::Vector{Vector{Int}}
    const prazo::Float64
    const coloracao::Vector{Int}
    const contagem::Matrix{Int}          # contagem[j, w] = nº de vizinhos de w com cor j
    const coloridos::Vector{Int}
    const proibido_k::BitVector
    const zobrist::Matrix{UInt128}       # zobrist[j, v]: chave de "v tem a cor j"
    const falhas::Set{UInt128}
    const max_transposicoes::Int
    const min_subarvore::Int
    hash::UInt128
    # áreas de trabalho da poda de Hall
    const par_do_vizinho::Vector{Int}    # par_do_vizinho[i] = cor emparelhada ao i-ésimo vizinho
    const visitado::BitVector
    const faltantes::Vector{Int}
    # estatísticas
    nos::Int
    podas_tabela::Int
    podas_hall::Int
end

"""
    EstatisticasV5

Contadores acumulados de uma execução do `bbg_v5`: nós da busca, ramos podados pela
tabela de transposição e pela poda de Hall, e estados guardados na tabela.
"""
mutable struct EstatisticasV5
    nos::Int
    podas_tabela::Int
    podas_hall::Int
    estados_guardados::Int
end
EstatisticasV5() = EstatisticasV5(0, 0, 0, 0)

"""
    bbg_v5(g; limite_superior=nothing, tempo_limite=900.0, max_transposicoes=2_000_000,
           min_subarvore=32, transposicao=true, hall=false, estatisticas=nothing) -> Resultado

Número de Grundy exato pela busca de k-átomos com limite por aresta (v4), tabela de
transposição e, opcionalmente, poda de Hall. `transposicao` e `hall` ligam/desligam
cada melhoria (para medir o efeito de cada uma); a poda de Hall vem desligada por
padrão porque, nas medições, custou mais do que economizou. Passe `estatisticas = EstatisticasV5()`
para receber os contadores da busca.
"""
function bbg_v5(g::SimpleGraph; limite_superior::Union{Nothing,Integer}=nothing,
                tempo_limite::Real=900.0, max_transposicoes::Integer=2_000_000,
                min_subarvore::Integer=32,
                transposicao::Bool=true, hall::Bool=false,
                estatisticas::Union{Nothing,EstatisticasV5}=nothing,
                rng::AbstractRNG=MersenneTwister(0x5eed))
    inicio = time()
    nv(g) == 0 && return Resultado(0, true, 0.0, Int[])
    b, entrada = limites_por_aresta(g)
    b = min.(b, limites_por_vertice(g))
    est = estatisticas === nothing ? EstatisticasV5() : estatisticas
    existe = (g, k, b, prazo; entrada) ->
        _existe_coloracao_parcial_v5(g, k, b, prazo, entrada, rng, est;
                                     max_transposicoes=transposicao ? Int(max_transposicoes) : 0,
                                     min_subarvore=Int(min_subarvore), hall)
    return _grundy_por_testemunhas(g, b, entrada, limite_superior, tempo_limite, inicio; existe)
end

function _existe_coloracao_parcial_v5(g::SimpleGraph, k::Int, b::Vector{Int}, prazo::Float64,
                                      entrada::Vector{Vector{Int}}, rng::AbstractRNG,
                                      est::EstatisticasV5; max_transposicoes::Int,
                                      min_subarvore::Int, hall::Bool)
    n = nv(g)
    Δg = n == 0 ? 0 : Δ(g)
    estado = EstadoV5(g, k, b, entrada, prazo, fill(SEM_COR, n), zeros(Int, k, n), Int[], falses(n),
                      rand(rng, UInt128, k, n), Set{UInt128}(), max_transposicoes, min_subarvore, UInt128(0),
                      zeros(Int, Δg), falses(k), Int[], 0, 0, 0)
    raizes = sort!(filter(v -> b[v] >= k, collect(vertices(g))); by=v -> (b[v], degree(g, v)), rev=true)
    try
        for r in raizes
            _colorir_v5!(estado, r, k)
            _buscar_v5!(estado, hall) && return copy(estado.coloracao)
            _descolorir_v5!(estado, r, k)
            estado.proibido_k[r] = true
        end
        return nothing
    finally
        est.nos += estado.nos
        est.podas_tabela += estado.podas_tabela
        est.podas_hall += estado.podas_hall
        est.estados_guardados += length(estado.falhas)
    end
end

@inline function _pode_receber(estado::EstadoV5, w::Int, i::Int, u::Int, j::Int)
    estado.coloracao[u] == SEM_COR && estado.entrada[w][i] >= j && estado.contagem[j, u] == 0 &&
        !(j == estado.k && estado.proibido_k[u])
end

function _buscar_v5!(estado::EstadoV5, hall::Bool)
    estado.nos += 1
    (estado.nos & 0x3ff) == 0 && _verificar_prazo(estado.prazo)

    if estado.max_transposicoes > 0 && estado.hash in estado.falhas
        estado.podas_tabela += 1
        return false
    end
    nos_antes = estado.nos
    resultado = _expandir_v5!(estado, hall)
    if !resultado && estado.nos - nos_antes >= estado.min_subarvore &&
       length(estado.falhas) < estado.max_transposicoes
        push!(estado.falhas, estado.hash)
    end
    return resultado
end

function _expandir_v5!(estado::EstadoV5, hall::Bool)
    g, contagem = estado.g, estado.contagem

    # Demanda mais restrita (MRV), como na v3/v4.
    melhor_w, melhor_j, melhor_qtd = 0, 0, typemax(Int)
    for w in estado.coloridos
        c = estado.coloracao[w]
        faltam = 0
        menor_qtd_w = typemax(Int)
        for j in 1:c-1
            contagem[j, w] > 0 && continue
            faltam += 1
            qtd = 0
            for (i, u) in enumerate(neighbors(g, w))
                _pode_receber(estado, w, i, u, j) && (qtd += 1)
            end
            qtd == 0 && return false
            menor_qtd_w = min(menor_qtd_w, qtd)
            qtd < melhor_qtd && ((melhor_w, melhor_j, melhor_qtd) = (w, j, qtd))
        end
        # Se toda cor faltante tem ≥ `faltam` candidatos, o emparelhamento existe
        # (Hall: qualquer conjunto S de cores tem ≥ faltam ≥ |S| candidatos) e o
        # algoritmo de Kuhn não precisa rodar.
        if hall
            if faltam >= 2 && menor_qtd_w < faltam && !_hall_satisfeito(estado, w)
                estado.podas_hall += 1
                return false
            end
        elseif faltam > 0
            livres = count(u -> estado.coloracao[u] == SEM_COR, neighbors(g, w))
            faltam > livres && return false
        end
    end
    melhor_w == 0 && return true

    candidatos = [u for (i, u) in enumerate(neighbors(g, melhor_w)) if _pode_receber(estado, melhor_w, i, u, melhor_j)]
    sort!(candidatos; by=u -> (estado.b[u], degree(g, u)), rev=true)
    for u in candidatos
        _colorir_v5!(estado, u, melhor_j)
        _buscar_v5!(estado, hall) && return true
        _descolorir_v5!(estado, u, melhor_j)
    end
    return false
end

# Existe emparelhamento que dá a cada cor faltante de w um vizinho distinto sem cor
# que pode recebê-la? (Kuhn: um caminho aumentante por cor.)
function _hall_satisfeito(estado::EstadoV5, w::Int)
    g, contagem = estado.g, estado.contagem
    faltantes = estado.faltantes
    empty!(faltantes)
    for j in 1:estado.coloracao[w]-1
        contagem[j, w] == 0 && push!(faltantes, j)
    end
    d = degree(g, w)
    par = estado.par_do_vizinho
    fill!(view(par, 1:d), 0)
    for j in faltantes
        fill!(estado.visitado, false)
        _aumentar!(estado, w, j) || return false
    end
    return true
end

function _aumentar!(estado::EstadoV5, w::Int, j::Int)
    estado.visitado[j] && return false
    estado.visitado[j] = true
    par = estado.par_do_vizinho
    for (i, u) in enumerate(neighbors(estado.g, w))
        _pode_receber(estado, w, i, u, j) || continue
        if par[i] == 0 || _aumentar!(estado, w, par[i])
            par[i] = j
            return true
        end
    end
    return false
end

function _colorir_v5!(estado::EstadoV5, v::Int, cor::Int)
    estado.coloracao[v] = cor
    push!(estado.coloridos, v)
    estado.hash ⊻= estado.zobrist[cor, v]
    for u in neighbors(estado.g, v)
        estado.contagem[cor, u] += 1
    end
end

function _descolorir_v5!(estado::EstadoV5, v::Int, cor::Int)
    for u in neighbors(estado.g, v)
        estado.contagem[cor, u] -= 1
    end
    estado.hash ⊻= estado.zobrist[cor, v]
    pop!(estado.coloridos)
    estado.coloracao[v] = SEM_COR
end
