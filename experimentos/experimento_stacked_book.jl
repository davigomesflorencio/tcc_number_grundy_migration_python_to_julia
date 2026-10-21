# Testes computacionais em grafos da classe Stacked Book (K_{1,m} □ P_n).
#
# Uso:
#   julia --project=. experimentos/experimento_stacked_book.jl ALGORITMO [TEMPO_LIMITE] [SAIDA]
#
# Exemplos:
#   julia --project=. experimentos/experimento_stacked_book.jl PLI_6
#   julia --project=. experimentos/experimento_stacked_book.jl BB_STF 900 resultados/classe/BB_STF.csv
#
# ALGORITMO ∈ PDG, BB_DELTA, BB_STF, BB_6, PLI_N, PLI_DELTA, PLI_STF, PLI_6

include(joinpath(@__DIR__, "comum.jl"))

const NS = 2:9
const MS = 3:9

function main(args)
    nome = argumento(args, 1, "PLI_6")
    tempo_limite = parse(Float64, argumento(args, 2, "900"))
    saida = argumento(args, 3, joinpath(RAIZ, "resultados", "classe", "$nome.csv"))
    preparar_csv(saida, ["n", "m", "vertices", "arestas", "delta", "stairfactor", "densidade", "grundy", "tempo"])
    aquecer(nome)

    for n in NS, m in MS
        g = ler_gml(joinpath(RAIZ, "instancias", "classe", "classe_graph_$(n)_$(m)"))
        V, E, delta, stf, dens = metricas(g)
        r, tempo = executar(nome, g, tempo_limite)
        anexar_linha(saida, (n, m, V, E, delta, stf, dens, r.grundy, tempo))
        @printf("%-9s n=%d m=%d |V|=%-3d grundy=%-3d tempo=%.4fs%s\n",
                nome, n, m, V, r.grundy, tempo, r.otimo ? "" : "  (limite de tempo)")
    end
    println("Resultados em $saida")
end

main(ARGS)
