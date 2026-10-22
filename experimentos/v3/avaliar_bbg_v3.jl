# Executa o BBG v3 nas instâncias do TCC e compara com o melhor valor que o TCC
# obteve (maior `grundy` entre todos os CSVs de `dados/` ou `dados_random/`).
#
# Uso:
#   julia --project=. experimentos/v3/avaliar_bbg_v3.jl [TEMPO_LIMITE] [SAIDA]

include(joinpath(@__DIR__, "..", "comum.jl"))

# melhor valor do TCC por instância, e se algum algoritmo do TCC a resolveu (< 900 s)
function melhores_tcc(diretorio, chave)
    melhor = Dict{String,Tuple{Int,Bool}}()
    for arquivo in filter(endswith(".csv"), readdir(joinpath(RAIZ, diretorio); join=true))
        linhas = readlines(arquivo)
        cab = split(linhas[1], ",")
        contador = Dict{String,Int}()
        for l in linhas[2:end]
            isempty(l) && continue
            d = Dict(zip(cab, split(l, ",")))
            id = chave(d, contador)
            valor = round(Int, parse(Float64, d["grundy"]))
            resolvido = parse(Float64, d["tempo"]) < 900
            v0, r0 = get(melhor, id, (0, false))
            melhor[id] = (max(v0, valor), r0 || resolvido)
        end
    end
    return melhor
end

chave_classe(d, _) = "classe_$(d["n"])_$(d["m"])"
function chave_random(d, contador)
    base = "$(d["vertices"])_$(d["prob"])"
    i = get(contador, base, 0)
    contador[base] = i + 1
    return "random_$(i)_$(base)"
end

function main(args)
    tempo_limite = parse(Float64, argumento(args, 1, "30"))
    saida = argumento(args, 2, joinpath(RAIZ, "resultados", "v3", "bbg_v3.csv"))
    preparar_csv(saida, ["instancia", "vertices", "arestas", "stairfactor", "ub_v3",
                         "grundy_tcc", "tcc_resolveu", "grundy_v3", "otimo_v3", "tempo_v3"])
    tcc = merge(melhores_tcc("dados", chave_classe), melhores_tcc("dados_random", chave_random))
    bbg_v3(grafo_stacked_book(2, 3))              # aquecimento do JIT

    instancias = Tuple{String,String}[]
    for n in 2:9, m in 3:9
        push!(instancias, ("classe_$(n)_$(m)", joinpath(RAIZ, "instancias", "classe", "classe_graph_$(n)_$(m)")))
    end
    for n in (10, 15, 20, 25), p in (0.1, 0.3, 0.5, 0.7, 0.9), i in 0:4
        push!(instancias, ("random_$(i)_$(n)_$(p)", joinpath(RAIZ, "instancias", "random", "random_graph_$(i)_$(n)_$(p).gml")))
    end

    resolvidas = 0; resolvidas_tcc = 0; melhorou = 0; pior = 0
    for (nome, caminho) in instancias
        g = ler_gml(caminho)
        gt, rt = get(tcc, nome, (0, false))
        r = bbg_v3(g; tempo_limite)
        anexar_linha(saida, (nome, nv(g), ne(g), stair_factor(g), limite_superior_vertices(g),
                             gt, rt, r.grundy, r.otimo, round(r.tempo; digits=4)))
        resolvidas += r.otimo; resolvidas_tcc += rt
        r.grundy > gt && (melhorou += 1)
        r.otimo && r.grundy < gt && (pior += 1)
        @printf("%-20s |V|=%-3d  TCC: %2d%s  v3: %2d %8.3fs%s\n", nome, nv(g),
                gt, rt ? " " : "*", r.grundy, r.tempo, r.otimo ? " " : "*")
    end
    println("\n* = não resolvida (TCC: 900 s em todos os algoritmos; v3: $(tempo_limite) s)")
    println("v3 resolveu $resolvidas/$(length(instancias)); o TCC resolveu $resolvidas_tcc/$(length(instancias))")
    println("v3 encontrou valor maior que o melhor do TCC em $melhorou instâncias")
    pior > 0 && println("ATENÇÃO: v3 provou ótimo menor que um valor do TCC em $pior instâncias")
    println("Resultados em $saida")
end

main(ARGS)
