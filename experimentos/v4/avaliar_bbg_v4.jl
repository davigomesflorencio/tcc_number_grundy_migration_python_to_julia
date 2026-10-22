# Executa o BBG v4 nas instâncias do TCC e compara com a v3 e com os valores
# concluídos do TCC.
#
# Uso:
#   julia --project=. experimentos/v4/avaliar_bbg_v4.jl [TEMPO_LIMITE] [SAIDA] [CSV_V3]
#
# CSV_V3 (padrão: resultados/v3/bbg_v3.csv) é a saída de `experimentos/v3/avaliar_bbg_v3.jl`;
# se existir, os resultados da v3 entram na comparação sem precisar rodá-la de novo.

include(joinpath(@__DIR__, "..", "comum.jl"))

function ler_resultados_v3(caminho)
    isfile(caminho) || return Dict{String,Dict{String,String}}()
    linhas = readlines(caminho)
    cab = split(linhas[1], ",")
    return Dict(String(split(l, ",")[1]) => Dict(String(k) => String(v) for (k, v) in zip(cab, split(l, ",")))
                for l in linhas[2:end] if !isempty(l))
end

function main(args)
    tempo_limite = parse(Float64, argumento(args, 1, "30"))
    saida = argumento(args, 2, joinpath(RAIZ, "resultados", "v4", "bbg_v4.csv"))
    v3 = ler_resultados_v3(argumento(args, 3, joinpath(RAIZ, "resultados", "v3", "bbg_v3.csv")))
    preparar_csv(saida, ["instancia", "vertices", "arestas", "ub_v3", "ub_v4", "grundy_tcc", "tcc_resolveu",
                         "grundy_v3", "otimo_v3", "tempo_v3", "grundy_v4", "otimo_v4", "tempo_v4"])
    tcc = valores_tcc()
    bbg_v4(grafo_stacked_book(2, 3))              # aquecimento do JIT

    resolvidas = Dict("v3" => 0, "v4" => 0, "tcc" => 0)
    so_v4 = 0; so_v3 = 0; ub_igual = 0; conflitos = 0
    for (nome, caminho) in instancias_tcc()
        g = ler_gml(caminho)
        t = get(tcc, nome, (confiavel=0, estourado=0, resolveu=false))
        r = bbg_v4(g; tempo_limite)
        ub4 = limite_superior_arestas(g)
        a3 = get(v3, nome, Dict{String,String}())
        ok3 = get(a3, "otimo_v3", "") == "true"
        anexar_linha(saida, (nome, nv(g), ne(g), limite_superior_vertices(g), ub4,
                             t.resolveu ? t.confiavel : "", t.resolveu,
                             get(a3, "grundy_v3", ""), get(a3, "otimo_v3", ""), get(a3, "tempo_v3", ""),
                             r.grundy, r.otimo, round(r.tempo; digits=4)))
        resolvidas["v4"] += r.otimo; resolvidas["v3"] += ok3; resolvidas["tcc"] += t.resolveu
        r.otimo && !ok3 && (so_v4 += 1)
        ok3 && !r.otimo && (so_v3 += 1)
        r.otimo && r.grundy == ub4 && (ub_igual += 1)
        t.resolveu && r.otimo && r.grundy != t.confiavel && (conflitos += 1)
        @printf("%-20s |V|=%-3d UB v3/v4=%2d/%-2d  TCC: %3s  v3: %3s%s  v4: %2d %8.3fs%s\n", nome, nv(g),
                limite_superior_vertices(g), ub4, t.resolveu ? string(t.confiavel) : "-",
                get(a3, "grundy_v3", "?"), ok3 ? " " : "*", r.grundy, r.tempo, r.otimo ? " " : "*")
    end
    total = length(instancias_tcc())
    println("\nTCC: - = nenhum algoritmo do TCC terminou em 900 s;  v3/v4: * = limite de $(tempo_limite) s")
    println("resolvidas: v4 $(resolvidas["v4"])/$total | v3 $(resolvidas["v3"])/$total | TCC $(resolvidas["tcc"])/$total")
    println("só a v4 resolveu: $so_v4 | só a v3 resolveu: $so_v3")
    println("v4 resolvidas com Γ = UB por aresta (sem precisar refutar): $ub_igual")
    conflitos > 0 && println("ATENÇÃO: v4 provou ótimo diferente de um resultado concluído do TCC em $conflitos instâncias")
    println("Resultados em $saida")
end

main(ARGS)
