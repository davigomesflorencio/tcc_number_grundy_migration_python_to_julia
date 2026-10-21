# Compara o BBG original (TCC) com o BBG v2 nas instâncias do TCC.
#
# Uso:
#   julia --project=. experimentos/v2/comparar_bbg.jl [LIMITE] [TEMPO_LIMITE] [SAIDA]
#
# LIMITE ∈ DELTA, STF, 6 (padrão: STF). Grava um CSV com uma linha por instância e
# imprime um resumo (instâncias resolvidas, tempo total, valores divergentes).

include(joinpath(@__DIR__, "..", "comum.jl"))

const LIMITES = Dict(
    "DELTA" => g -> grau_maximo(g) + 1,
    "STF"   => stair_factor,
    "6"     => g -> 6,
)

function instancias()
    lista = Tuple{String,Any}[]
    for n in 2:9, m in 3:9
        push!(lista, ("classe_$(n)_$(m)", ler_gml(joinpath(RAIZ, "instancias", "classe", "classe_graph_$(n)_$(m)"))))
    end
    for n in (10, 15, 20, 25), p in (0.1, 0.3, 0.5, 0.7, 0.9), i in 0:4
        push!(lista, ("random_$(i)_$(n)_$(p)", ler_gml(joinpath(RAIZ, "instancias", "random", "random_graph_$(i)_$(n)_$(p).gml"))))
    end
    return lista
end

function main(args)
    nome_limite = argumento(args, 1, "STF")
    tempo_limite = parse(Float64, argumento(args, 2, "60"))
    saida = argumento(args, 3, joinpath(RAIZ, "resultados", "v2", "comparacao_bbg_$(nome_limite).csv"))
    limite = LIMITES[nome_limite]
    preparar_csv(saida, ["instancia", "vertices", "arestas", "ub",
                         "grundy_v1", "otimo_v1", "tempo_v1", "grundy_v2", "otimo_v2", "tempo_v2"])

    g0 = grafo_stacked_book(2, 3)                 # aquecimento do JIT
    bbg(g0, 4); bbg_v2(g0, 4)

    resolvidos = [0, 0]; tempo_total = [0.0, 0.0]; divergentes = 0
    for (nome, g) in instancias()
        ub = limite(g)
        r1 = bbg(g, ub; tempo_limite)
        r2 = bbg_v2(g, ub; tempo_limite)
        anexar_linha(saida, (nome, nv(g), ne(g), ub, r1.grundy, r1.otimo, round(r1.tempo; digits=4),
                             r2.grundy, r2.otimo, round(r2.tempo; digits=4)))
        resolvidos .+= (r1.otimo, r2.otimo)
        tempo_total .+= (r1.tempo, r2.tempo)
        r1.grundy != r2.grundy && (divergentes += 1)
        @printf("%-20s |V|=%-3d UB=%-3d  v1: %2d %8.3fs%s   v2: %2d %8.3fs%s\n", nome, nv(g), ub,
                r1.grundy, r1.tempo, r1.otimo ? " " : "*", r2.grundy, r2.tempo, r2.otimo ? " " : "*")
    end
    total = length(instancias())
    println("\n* = limite de tempo ($(tempo_limite)s)")
    @printf("v1: %d/%d resolvidas, tempo total %.1fs\n", resolvidos[1], total, tempo_total[1])
    @printf("v2: %d/%d resolvidas, tempo total %.1fs\n", resolvidos[2], total, tempo_total[2])
    println("instâncias com valores diferentes: $divergentes")
    println("Resultados em $saida")
end

main(ARGS)
