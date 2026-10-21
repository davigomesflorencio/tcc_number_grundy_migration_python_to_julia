# Testes computacionais em grafos aleatórios G(n, p).
#
# Uso:
#   julia --project=. experimentos/experimento_random.jl ALGORITMO [TEMPO_LIMITE] [SAIDA]
#
# Exemplos:
#   julia --project=. experimentos/experimento_random.jl PLI_DELTA
#   julia --project=. experimentos/experimento_random.jl BB_STF 900 resultados/random/BB_STF.csv
#
# ALGORITMO ∈ PDG, BB_DELTA, BB_STF, BB_6, PLI_N, PLI_DELTA, PLI_STF, PLI_6

include(joinpath(@__DIR__, "comum.jl"))

const VERTICES = [10, 15, 20, 25]
const PROBABILIDADES = [0.1, 0.3, 0.5, 0.7, 0.9]
const REPETICOES = 0:4

function main(args)
    nome = argumento(args, 1, "PLI_DELTA")
    tempo_limite = parse(Float64, argumento(args, 2, "900"))
    saida = argumento(args, 3, joinpath(RAIZ, "resultados", "random", "$nome.csv"))
    preparar_csv(saida, ["prob", "vertices", "arestas", "delta", "stairfactor", "densidade", "grundy", "tempo"])
    aquecer(nome)

    for n in VERTICES, p in PROBABILIDADES, i in REPETICOES
        g = ler_gml(joinpath(RAIZ, "instancias", "random", "random_graph_$(i)_$(n)_$(p).gml"))
        V, E, delta, stf, dens = metricas(g)
        r, tempo = executar(nome, g, tempo_limite)
        anexar_linha(saida, (p, V, E, delta, stf, dens, Float64(r.grundy), tempo))
        @printf("%-9s n=%-3d p=%.1f i=%d  grundy=%-3d tempo=%.4fs%s\n",
                nome, n, p, i, r.grundy, tempo, r.otimo ? "" : "  (limite de tempo)")
    end
    println("Resultados em $saida")
end

main(ARGS)
