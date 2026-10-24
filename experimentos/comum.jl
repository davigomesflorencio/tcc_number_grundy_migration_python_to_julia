# Utilitários compartilhados pelos scripts de experimentos computacionais.

using NumberGrundy
using Graphs
using Printf

const RAIZ = normpath(joinpath(@__DIR__, ".."))

"""
Algoritmos avaliados no TCC, identificados pelo nome usado nos arquivos CSV.
Cada entrada recebe `(g, tempo_limite)` e devolve um `Resultado`.
"""
const ALGORITMOS = Dict{String,Function}(
    "PDG"       => (g, t) -> pdg(g; tempo_limite=t),
    "BB_DELTA"  => (g, t) -> bbg(g, grau_maximo(g) + 1; tempo_limite=t),
    "BB_STF"    => (g, t) -> bbg(g, stair_factor(g); tempo_limite=t),
    "BB_6"      => (g, t) -> bbg(g, 6; tempo_limite=t),
    "PLI_N"     => (g, t) -> pli_grundy(g, :n; tempo_limite=t),
    "PLI_DELTA" => (g, t) -> pli_grundy(g, :delta; tempo_limite=t),
    "PLI_STF"   => (g, t) -> pli_grundy(g, :stf; tempo_limite=t),
    "PLI_6"     => (g, t) -> pli_grundy(g, :prod; tempo_limite=t),
    # v2
    "BB2_DELTA" => (g, t) -> bbg_v2(g, grau_maximo(g) + 1; tempo_limite=t),
    "BB2_STF"   => (g, t) -> bbg_v2(g, stair_factor(g); tempo_limite=t),
    "BB2_6"     => (g, t) -> bbg_v2(g, 6; tempo_limite=t),
    # v3
    "BB3"       => (g, t) -> bbg_v3(g; tempo_limite=t),
    # v4
    "BB4"       => (g, t) -> bbg_v4(g; tempo_limite=t),
    "PLI_ARESTA" => (g, t) -> pli_grundy(g, limite_superior_arestas(g); tempo_limite=t),
    # v5
    "BB5"       => (g, t) -> bbg_v5(g; tempo_limite=t),
    # formulação de representantes
    "PLI_REP"   => (g, t) -> pli_representantes(g, :stf; tempo_limite=t),
)

densidade(g) = nv(g) <= 1 ? 0.0 : 2ne(g) / (nv(g) * (nv(g) - 1))

"""
Métricas da instância na ordem das colunas dos CSVs:
vertices, arestas, delta, stairfactor, densidade.
"""
metricas(g) = (nv(g), ne(g), grau_maximo(g) + 1, stair_factor(g), round(densidade(g); digits=4))

function executar(nome::AbstractString, g, tempo_limite)
    haskey(ALGORITMOS, nome) || error("algoritmo '$nome' desconhecido. Opções: $(join(sort!(collect(keys(ALGORITMOS))), ", "))")
    r = ALGORITMOS[nome](g, tempo_limite)
    # como no TCC, uma execução que estoura o limite é registrada com o tempo limite
    tempo = r.otimo ? round(r.tempo; digits=4) : Float64(tempo_limite)
    return r, tempo
end

function anexar_linha(caminho, valores)
    open(caminho, "a") do io
        println(io, join(string.(valores), ","))
    end
end

function preparar_csv(caminho, cabecalho)
    mkpath(dirname(caminho))
    isfile(caminho) || open(io -> println(io, join(cabecalho, ",")), caminho, "w")
    return caminho
end

function argumento(args, i, padrao)
    i <= length(args) ? args[i] : padrao
end

# Aquece a compilação JIT para que o primeiro tempo medido não inclua compilação.
function aquecer(nome)
    g = grafo_stacked_book(2, 3)
    ALGORITMOS[nome](g, 60.0)
    return nothing
end

# ---------------------------------------------------------------------------
# Instâncias e valores de referência do TCC

"""
    instancias_tcc() -> Vector{Tuple{String,String}}

As 156 instâncias do TCC como pares (nome, caminho): 56 Stacked Book e 100 aleatórias.
"""
function instancias_tcc()
    instancias = Tuple{String,String}[]
    for n in 2:9, m in 3:9
        push!(instancias, ("classe_$(n)_$(m)", joinpath(RAIZ, "instancias", "classe", "classe_graph_$(n)_$(m)")))
    end
    for n in (10, 15, 20, 25), p in (0.1, 0.3, 0.5, 0.7, 0.9), i in 0:4
        push!(instancias, ("random_$(i)_$(n)_$(p)", joinpath(RAIZ, "instancias", "random", "random_graph_$(i)_$(n)_$(p).gml")))
    end
    return instancias
end

# Só são confiáveis os valores do TCC de execuções que terminaram antes de 900 s:
# nas que estouraram, o Python gravava `config.k`, que podia vir de outra instância
# (ver README, seção "Atenção: valores do TCC em execuções que estouraram 900 s").
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

valores_tcc() = merge(valores_tcc("dados", chave_classe), valores_tcc("dados_random", chave_random))
