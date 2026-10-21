using NumberGrundy
using Graphs
using Random
using Test

const RAIZ = normpath(joinpath(@__DIR__, ".."))

# Lê um CSV de resultados do TCC como vetor de Dict (coluna => texto).
function ler_csv(caminho)
    linhas = readlines(caminho)
    cab = split(linhas[1], ",")
    return [Dict(zip(cab, split(l, ","))) for l in linhas[2:end] if !isempty(l)]
end

# Instâncias aleatórias aparecem nos CSVs na ordem (n, p, i = 0..4).
function instancias_random(linhas)
    contador = Dict{Tuple{String,String},Int}()
    map(linhas) do l
        chave = (l["vertices"], l["prob"])
        i = get(contador, chave, 0)
        contador[chave] = i + 1
        ler_gml(joinpath(RAIZ, "instancias", "random", "random_graph_$(i)_$(l["vertices"])_$(l["prob"]).gml"))
    end
end

grundy(l) = round(Int, parse(Float64, l["grundy"]))
resolvido(l) = parse(Float64, l["tempo"]) < 900

@testset "NumberGrundy" begin
    @testset "Coloração gulosa" begin
        g = SimpleGraph(4)
        for (u, v) in [(1, 2), (1, 3), (2, 4), (3, 4), (2, 3)]
            add_edge!(g, u, v)
        end
        k, c = coloracao_gulosa(g, [1, 2, 3, 4])
        @test k == 3
        @test eh_coloracao_grundy(g, c)
        @test eh_coloracao_grundy(path_graph(3), [2, 1, 2])
        @test !eh_coloracao_grundy(path_graph(3), [1, 3, 2])   # vértice 3 (cor 2) sem vizinho de cor 1
        @test !eh_coloracao_grundy(path_graph(2), [1, 3])
        @test cor_minima(g, 4, [1, 2, 0, 0]) == 1
    end

    @testset "Cliques maximais" begin
        g = complete_graph(4)
        @test cliques_maximais(g) == [[1, 2, 3, 4]]
        p = path_graph(4)
        @test sort(cliques_maximais(p)) == [[1, 2], [2, 3], [3, 4]]
        # conjuntos independentes maximais do caminho P4: {1,3}, {1,4}, {2,4}
        adj = NumberGrundy.mascaras_adjacencia(p)
        encontrados = UInt64[]
        conjuntos_independentes_maximais(I -> push!(encontrados, I), adj, UInt64(0b1111))
        @test sort(encontrados) == sort(UInt64[0b0101, 0b1001, 0b1010])
    end

    @testset "Valores conhecidos" begin
        # Γ(P4) = 3, Γ(K_n) = n, Γ(K_{1,m}) = 2, Γ(C4) = 2
        for (g, esperado) in [(path_graph(4), 3), (complete_graph(5), 5), (star_graph(6), 2), (cycle_graph(4), 2)]
            @test pdg(g).grundy == esperado
            @test bbg(g, nv(g)).grundy == esperado
            @test pli_grundy(g, :n).grundy == esperado
        end
        @test pli_numero_cromatico(cycle_graph(5)).grundy == 3
        @test lb_deg(path_graph(4))[1] <= 3
    end

    @testset "Gerador Stacked Book = instâncias do TCC" begin
        for n in 2:9, m in 3:9
            @test grafo_stacked_book(n, m) == ler_gml(joinpath(RAIZ, "instancias", "classe", "classe_graph_$(n)_$(m)"))
        end
    end

    @testset "Escrita/leitura GML" begin
        g = grafo_aleatorio(12, 0.4; rng=MersenneTwister(1))
        caminho = escrever_gml(tempname(), g)
        @test ler_gml(caminho) == g
    end

    @testset "Limites (Δ+1, Stair Factor) = TCC" begin
        for l in ler_csv(joinpath(RAIZ, "dados", "PDG.csv"))
            g = ler_gml(joinpath(RAIZ, "instancias", "classe", "classe_graph_$(l["n"])_$(l["m"])"))
            @test grau_maximo(g) + 1 == parse(Int, l["delta"])
            @test stair_factor(g) == parse(Int, l["stairfactor"])
        end
        linhas = ler_csv(joinpath(RAIZ, "dados_random", "PLI_STF.csv"))
        for (l, g) in zip(linhas, instancias_random(linhas))
            @test ne(g) == parse(Int, l["arestas"])
            @test stair_factor(g) == parse(Int, l["stairfactor"])
        end
    end

    @testset "Número de Grundy = TCC (instâncias pequenas)" begin
        # Stacked Book com até 24 vértices: PDG, BBG e PLI
        for l in ler_csv(joinpath(RAIZ, "dados", "PDG.csv"))
            resolvido(l) || continue
            g = ler_gml(joinpath(RAIZ, "instancias", "classe", "classe_graph_$(l["n"])_$(l["m"])"))
            nv(g) <= 24 || continue
            r = pdg(g)
            @test r.grundy == grundy(l)
            @test eh_coloracao_grundy(g, r.coloracao) && maximum(r.coloracao) == r.grundy
            @test bbg(g, stair_factor(g)).grundy == grundy(l)
            @test pli_grundy(g, :stf).grundy == grundy(l)
        end
        # Aleatórias com 10 vértices
        linhas = filter(l -> l["vertices"] == "10", ler_csv(joinpath(RAIZ, "dados_random", "PLI_STF.csv")))
        for (l, g) in zip(linhas, instancias_random(linhas))
            @test pli_grundy(g, :delta).grundy == grundy(l)
            r = bbg(g, grau_maximo(g) + 1)
            @test r.grundy == grundy(l)
            @test eh_coloracao_grundy(g, r.coloracao)
            @test pdg(g).grundy == grundy(l)
        end
    end

    @testset "Limite de tempo devolve melhor valor parcial" begin
        g = ler_gml(joinpath(RAIZ, "instancias", "random", "random_graph_0_30_0.5.gml"))
        r = bbg(g, nv(g); tempo_limite=0.5)
        @test !r.otimo
        @test r.grundy == 0 || eh_coloracao_grundy(g, r.coloracao)
    end
end
