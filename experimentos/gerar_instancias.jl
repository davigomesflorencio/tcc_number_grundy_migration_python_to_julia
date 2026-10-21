# Gera novas instâncias (aleatórias e Stacked Book) no formato GML.
#
# Uso:
#   julia --project=. experimentos/gerar_instancias.jl [DIRETORIO_SAIDA] [SEMENTE]
#
# Por padrão grava em `instancias_geradas/` para não sobrescrever as instâncias
# originais do TCC em `instancias/`.

include(joinpath(@__DIR__, "comum.jl"))
using Random

function main(args)
    destino = argumento(args, 1, joinpath(RAIZ, "instancias_geradas"))
    rng = MersenneTwister(parse(Int, argumento(args, 2, "2022")))
    mkpath(joinpath(destino, "random"))
    mkpath(joinpath(destino, "classe"))

    for i in 0:9, n in (10, 15, 20, 25, 30), p in (0.1, 0.3, 0.5, 0.7, 0.9)
        escrever_gml(joinpath(destino, "random", "random_graph_$(i)_$(n)_$(p).gml"), grafo_aleatorio(n, p; rng))
    end
    for n in 2:9, m in 3:9
        escrever_gml(joinpath(destino, "classe", "classe_graph_$(n)_$(m)"), grafo_stacked_book(n, m))
    end
    println("Instâncias geradas em $destino")
end

main(ARGS)
