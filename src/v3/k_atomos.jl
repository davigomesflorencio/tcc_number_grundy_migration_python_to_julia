# BBG v3 — número de Grundy por busca de testemunhas (k-átomos).
#
# Coloração parcial de Grundy: coloração própria de um subconjunto S ⊆ V em que todo
# vértice de S com cor c tem, em S, vizinhos com todas as cores 1..c-1.
#
# Fato central: se existe uma coloração parcial de Grundy que usa a cor k, então
# Γ(G) ≥ k. Basta aplicar a coloração gulosa com os vértices de S primeiro (em
# ordem crescente de cor) e depois os demais: cada vértice de S recebe exatamente
# a sua cor. Reciprocamente, toda coloração de Grundy com k cores é uma coloração
# parcial de Grundy. Logo
#
#   Γ(G) ≥ k  ⇔  existe coloração parcial de Grundy com algum vértice de cor k.
#
# A busca não colore o grafo inteiro: escolhe uma raiz r com cor k e, a cada passo,
# uma "demanda" (vértice w de cor c sem vizinho de cor j < c) e ramifica sobre qual
# vizinho de w recebe a cor j. Cada vértice de cor c gera no máximo c-1 demandas,
# então a testemunha tem no máximo T(k) = 1 + Σ_{j<k} T(j) = 2^(k-1) vértices
# (o tamanho de um k-átomo, Zaker 2006), independentemente de n.
#
# Complexidade do teste "Γ ≥ k?": O(n · k · W · Δ^W), com W = min(n, 2^(k-1)):
# n raízes, profundidade ≤ W, ramificação ≤ Δ e O(W · k · Δ) por nó. Polinomial
# para k fixo. O algoritmo testa k = LB+1, LB+2, ... até a primeira resposta
# negativa (k = Γ+1), que domina o total: O(n · Γ · W · Δ^W) com W = min(n, 2^Γ).

mutable struct EstadoKAtomo
    const g::SimpleGraph{Int}
    const k::Int
    const b::Vector{Int}                 # limite da cor de cada vértice
    const entrada::Vector{Vector{Int}}   # entrada[w][i]: limite da cor do i-ésimo vizinho de w
                                         # quando ele é testemunha de w (v3: b do vizinho)
    const prazo::Float64
    const coloracao::Vector{Int}
    const contagem::Matrix{Int}          # contagem[j, w] = nº de vizinhos de w com cor j
    const coloridos::Vector{Int}
    const proibido_k::BitVector          # raízes já descartadas não recebem a cor k
    nos::Int
end

"""
    bbg_v3(g; limite_superior=nothing, tempo_limite=900.0) -> Resultado

Número de Grundy exato pela busca de colorações parciais de Grundy (k-átomos).
O limite superior padrão é `limite_superior_vertices(g)`, que é sempre válido; um
`limite_superior` informado é combinado com ele pelo mínimo.
"""
function bbg_v3(g::SimpleGraph; limite_superior::Union{Nothing,Integer}=nothing,
                tempo_limite::Real=900.0)
    inicio = time()
    nv(g) == 0 && return Resultado(0, true, 0.0, Int[])
    b = limites_por_vertice(g)
    return _grundy_por_testemunhas(g, b, entrada_por_vertice(g, b), limite_superior, tempo_limite, inicio)
end

"""
    entrada_por_vertice(g, b) -> Vector{Vector{Int}}

Limites de "entrada" sem informação de aresta: o i-ésimo vizinho `u` de `w` tem
limite `b[u]`.
"""
entrada_por_vertice(g::SimpleGraph, b::Vector{Int}) = [b[neighbors(g, w)] for w in vertices(g)]

# Laço externo comum à v3 e à v4: testa Γ ≥ LB+1, LB+2, ... até a primeira negativa.
function _grundy_por_testemunhas(g::SimpleGraph, b::Vector{Int}, entrada::Vector{Vector{Int}},
                                 limite_superior, tempo_limite::Real, inicio::Float64=time())
    prazo = inicio + tempo_limite
    ub = maximum(b)
    limite_superior === nothing || (ub = min(ub, Int(limite_superior)))

    lb, melhor = lb_deg(g)
    k_dec, c_dec = coloracao_gulosa(g, decomposicao_vertices(g)[1])
    k_dec > lb && ((lb, melhor) = (k_dec, c_dec))

    try
        while lb < ub
            parcial = _existe_coloracao_parcial(g, lb + 1, b, prazo; entrada)
            parcial === nothing && break                      # Γ(G) < lb + 1
            k, coloracao = _completar_gulosamente(g, parcial)
            (lb, melhor) = (k, coloracao)
        end
    catch e
        e isa TempoEsgotado || rethrow()
        return Resultado(lb, false, time() - inicio, melhor)
    end
    return Resultado(lb, true, time() - inicio, melhor)
end

# Estende uma coloração parcial de Grundy a uma coloração de Grundy de G.
function _completar_gulosamente(g::SimpleGraph, parcial::Vector{Int})
    coloridos = sort!(filter(v -> parcial[v] != SEM_COR, collect(vertices(g))); by=v -> parcial[v])
    resto = filter(v -> parcial[v] == SEM_COR, collect(vertices(g)))
    return coloracao_gulosa(g, vcat(coloridos, resto))
end

"""
    _existe_coloracao_parcial(g, k, b, prazo; entrada) -> Union{Vector{Int}, Nothing}

Coloração parcial de Grundy em que algum vértice tem a cor `k`, ou `nothing`.
"""
function _existe_coloracao_parcial(g::SimpleGraph, k::Int, b::Vector{Int}, prazo::Float64;
                                   entrada::Vector{Vector{Int}}=entrada_por_vertice(g, b))
    n = nv(g)
    estado = EstadoKAtomo(g, k, b, entrada, prazo, fill(SEM_COR, n), zeros(Int, k, n), Int[], falses(n), 0)
    raizes = sort!(filter(v -> b[v] >= k, collect(vertices(g))); by=v -> (b[v], degree(g, v)), rev=true)
    for r in raizes
        _colorir!(estado, r, k)
        _buscar!(estado) && return copy(estado.coloracao)
        _descolorir!(estado, r, k)
        estado.proibido_k[r] = true      # soluções com r de cor k já foram descartadas
    end
    return nothing
end

# `u` (o i-ésimo vizinho de w) pode receber a cor j como testemunha de w?
@inline function _pode_receber(estado::EstadoKAtomo, w::Int, i::Int, u::Int, j::Int)
    estado.coloracao[u] == SEM_COR && estado.entrada[w][i] >= j && estado.contagem[j, u] == 0 &&
        !(j == estado.k && estado.proibido_k[u])
end

function _buscar!(estado::EstadoKAtomo)
    estado.nos += 1
    (estado.nos & 0x3ff) == 0 && _verificar_prazo(estado.prazo)
    g, contagem = estado.g, estado.contagem

    # Escolhe a demanda (w, j) com menos candidatos (MRV). Poda se alguma não tem
    # candidato ou se w tem mais cores faltando do que vizinhos sem cor.
    melhor_w, melhor_j, melhor_qtd = 0, 0, typemax(Int)
    for w in estado.coloridos
        c = estado.coloracao[w]
        faltam = 0
        for j in 1:c-1
            contagem[j, w] > 0 && continue
            faltam += 1
            qtd = 0
            for (i, u) in enumerate(neighbors(g, w))
                _pode_receber(estado, w, i, u, j) && (qtd += 1)
            end
            qtd == 0 && return false
            qtd < melhor_qtd && ((melhor_w, melhor_j, melhor_qtd) = (w, j, qtd))
        end
        if faltam > 0
            livres = count(u -> estado.coloracao[u] == SEM_COR, neighbors(g, w))
            faltam > livres && return false
        end
    end
    melhor_w == 0 && return true             # nenhuma demanda: coloração parcial de Grundy

    candidatos = [u for (i, u) in enumerate(neighbors(g, melhor_w)) if _pode_receber(estado, melhor_w, i, u, melhor_j)]
    sort!(candidatos; by=u -> (estado.b[u], degree(g, u)), rev=true)
    for u in candidatos
        _colorir!(estado, u, melhor_j)
        _buscar!(estado) && return true
        _descolorir!(estado, u, melhor_j)
    end
    return false
end

function _colorir!(estado::EstadoKAtomo, v::Int, cor::Int)
    estado.coloracao[v] = cor
    push!(estado.coloridos, v)
    for u in neighbors(estado.g, v)
        estado.contagem[cor, u] += 1
    end
end

function _descolorir!(estado::EstadoKAtomo, v::Int, cor::Int)
    for u in neighbors(estado.g, v)
        estado.contagem[cor, u] -= 1
    end
    pop!(estado.coloridos)
    estado.coloracao[v] = SEM_COR
end
