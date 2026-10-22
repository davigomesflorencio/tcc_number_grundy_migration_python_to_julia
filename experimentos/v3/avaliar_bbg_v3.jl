# Executa o BBG v3 nas instâncias do TCC e compara com os valores do TCC.
#
# Só são considerados confiáveis os valores do TCC de execuções que terminaram
# antes do limite de 900 s. Nas execuções que estouraram o tempo, o Python gravava
# `config.k`, uma variável global que só era zerada após uma execução bem-sucedida
# (e que o `model_grundy` nem atualizava); por isso esses valores podem ter vindo
# de outra instância. Ex.: BB_DELTA = 9 em 17 grafos aleatórios seguidos com n = 15,
# cujo Γ real é 6, 7 ou 8. Esses valores aparecem apenas na coluna informativa
# `grundy_tcc_estourado`.
#
# Uso:
#   julia --project=. experimentos/v3/avaliar_bbg_v3.jl [TEMPO_LIMITE] [SAIDA]

include(joinpath(@__DIR__, "..", "comum.jl"))

const LIMITE_TCC = 900.0

# Por instância: maior valor do TCC entre execuções concluídas (`confiavel`, 0 se
# nenhuma terminou) e maior valor entre execuções que estouraram o tempo (`estourado`).
function valores_tcc(diretorio, chave)
    valores = Dict{String,NamedTuple{(:confiavel, :estourado, :resolveu),Tuple{Int,Int,Bool}}}()
    for arquivo in filter(endswith(".csv"), readdir(joinpath(RAIZ, diretorio); join=true))
        linhas = readlines(arquivo)
        cab = split(linhas[1], ",")
        contador = Dict{String,Int}()
        for l in linhas[2:end]
            isempty(l) && continue
            d = Dict(zip(cab, split(l, ",")))
            id = chave(d, contador)
            valor = round(Int, parse(Float64, d["grundy"]))
            terminou = parse(Float64, d["tempo"]) < LIMITE_TCC
            a = get(valores, id, (confiavel=0, estourado=0, resolveu=false))
            valores[id] = terminou ?
                (confiavel=max(a.confiavel, valor), estourado=a.estourado, resolveu=true) :
                (confiavel=a.confiavel, estourado=max(a.estourado, valor), resolveu=a.resolveu)
        end
    end
    return valores
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
                         "grundy_tcc", "tcc_resolveu", "grundy_tcc_estourado", "grundy_v3", "otimo_v3", "tempo_v3"])
    tcc = merge(valores_tcc("dados", chave_classe), valores_tcc("dados_random", chave_random))
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
        t = get(tcc, nome, (confiavel=0, estourado=0, resolveu=false))
        gt, rt = t.confiavel, t.resolveu
        r = bbg_v3(g; tempo_limite)
        anexar_linha(saida, (nome, nv(g), ne(g), stair_factor(g), limite_superior_vertices(g),
                             rt ? gt : "", rt, t.estourado, r.grundy, r.otimo, round(r.tempo; digits=4)))
        resolvidas += r.otimo; resolvidas_tcc += rt
        rt && r.grundy > gt && (melhorou += 1)
        rt && r.otimo && r.grundy < gt && (pior += 1)
        @printf("%-20s |V|=%-3d  TCC: %3s  v3: %2d %8.3fs%s\n", nome, nv(g),
                rt ? string(gt) : "-", r.grundy, r.tempo, r.otimo ? " " : "*")
    end
    println("\nTCC: - = nenhum algoritmo do TCC terminou em 900 s;  v3: * = limite de $(tempo_limite) s")
    println("v3 resolveu $resolvidas/$(length(instancias)); o TCC resolveu $resolvidas_tcc/$(length(instancias))")
    println("v3 com valor maior que um resultado concluído do TCC: $melhorou instâncias")
    pior > 0 && println("ATENÇÃO: v3 provou ótimo menor que um resultado concluído do TCC em $pior instâncias")
    println("Resultados em $saida")
end

main(ARGS)
