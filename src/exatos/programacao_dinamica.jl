# PDG — Algoritmo de Programação Dinâmica para o número de Grundy
# (baseado no algoritmo de Lawler para coloração).
#
# Para S ⊆ V:  Γ(∅) = 0  e  Γ(S) = max { Γ(S \ I) + 1 : I conjunto independente
# maximal de G[S] }.  Os conjuntos são representados por bitmask e percorridos em
# ordem crescente, de forma que todo subconjunto de S já foi calculado antes de S
# (mesma ordem do `powerset` do Python). (Python: `lawler`.)

"""
    pdg(g; tempo_limite=900.0, n_max=30) -> Resultado

Número de Grundy exato via Programação Dinâmica sobre subconjuntos. Memória
O(2^n) bytes, por isso `nv(g)` é limitado a `n_max`. Se o tempo esgotar, devolve o
maior Γ(S) já calculado (limite inferior válido, pois Γ é monótono em subgrafos
induzidos).
"""
function pdg(g::SimpleGraph; tempo_limite::Real=900.0, n_max::Integer=30)
    n = nv(g)
    n <= n_max || throw(ArgumentError("pdg: n = $n excede n_max = $n_max (memória 2^n)"))
    inicio = time()
    prazo = inicio + tempo_limite
    adj = mascaras_adjacencia(g)
    total = UInt64(1) << n
    gama = zeros(UInt8, total)        # gama[S+1] = Γ(G[S])
    melhor = 0

    try
        for S in UInt64(1):total-1
            (S & 0x3ff) == 0 && _verificar_prazo(prazo)
            valor = 0
            conjuntos_independentes_maximais(adj, S) do I
                valor = max(valor, gama[(S & ~I)+1] + 1)
            end
            gama[S+1] = valor
            melhor = max(melhor, valor)
        end
    catch e
        e isa TempoEsgotado || rethrow()
        return Resultado(melhor, false, time() - inicio, Int[])
    end

    coloracao = _reconstruir_coloracao(adj, gama, total - 1, n)
    return Resultado(Int(gama[end]), true, time() - inicio, coloracao)
end

# Refaz o caminho ótimo da PD: a cor 1 é o conjunto I escolhido em S, a cor 2 é o
# conjunto escolhido em S \ I, e assim por diante.
function _reconstruir_coloracao(adj, gama, S::UInt64, n)
    coloracao = zeros(Int, n)
    cor = 1
    while S != 0
        alvo = gama[S+1] - 1
        escolhido = UInt64(0)
        conjuntos_independentes_maximais(adj, S) do I
            escolhido == 0 && gama[(S & ~I)+1] == alvo && (escolhido = I)
        end
        for v in 1:n
            (escolhido >> (v - 1)) & 1 == 1 && (coloracao[v] = cor)
        end
        S &= ~escolhido
        cor += 1
    end
    return coloracao
end
