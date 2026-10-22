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
