# Leitura e escrita de grafos no formato GML gerado pelo NetworkX (`nx.write_gml`).
#
# As instâncias do TCC possuem nós com `id` e `label` inteiros (0..n-1) e arestas
# `source`/`target` referenciando os `id`s. Equivalente a
# `nx.read_gml(caminho, destringizer=int)`: o rótulo `k` vira o vértice `k+1`.

"""
    ler_gml(caminho) -> SimpleGraph

Lê uma instância GML do NetworkX. O vértice Julia `v` corresponde ao rótulo `v-1`.
"""
function ler_gml(caminho::AbstractString)
    texto = read(caminho, String)
    rotulo_por_id = Dict{Int,Int}()
    arestas = Tuple{Int,Int}[]

    for bloco in eachmatch(r"node\s*\[(.*?)\]"s, texto)
        id = parse(Int, match(r"\bid\s+(-?\d+)", bloco[1])[1])
        m = match(r"\blabel\s+\"?(-?\d+)\"?", bloco[1])
        rotulo_por_id[id] = m === nothing ? id : parse(Int, m[1])
    end
    for bloco in eachmatch(r"edge\s*\[(.*?)\]"s, texto)
        s = parse(Int, match(r"\bsource\s+(-?\d+)", bloco[1])[1])
        t = parse(Int, match(r"\btarget\s+(-?\d+)", bloco[1])[1])
        push!(arestas, (s, t))
    end

    rotulos = sort!(collect(values(rotulo_por_id)))
    rotulos == collect(0:length(rotulos)-1) ||
        error("ler_gml: rótulos devem ser os inteiros 0..n-1 ($caminho)")

    g = SimpleGraph(length(rotulos))
    for (s, t) in arestas
        add_edge!(g, rotulo_por_id[s] + 1, rotulo_por_id[t] + 1)
    end
    return g
end

"""
    escrever_gml(caminho, g)

Escreve `g` no mesmo formato GML produzido pelo NetworkX (rótulos 0..n-1).
"""
function escrever_gml(caminho::AbstractString, g::SimpleGraph)
    open(caminho, "w") do io
        println(io, "graph [")
        for v in vertices(g)
            println(io, "  node [\n    id $(v-1)\n    label \"$(v-1)\"\n  ]")
        end
        for e in edges(g)
            println(io, "  edge [\n    source $(src(e)-1)\n    target $(dst(e)-1)\n  ]")
        end
        println(io, "]")
    end
    return caminho
end
