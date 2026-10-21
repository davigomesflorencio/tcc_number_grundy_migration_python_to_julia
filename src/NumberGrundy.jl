"""
    NumberGrundy

Algoritmos para o Problema do Número de Grundy (número cromático guloso), migrados
do TCC original em Python (NetworkX + OR-Tools) para Julia (Graphs.jl + JuMP).

Organização:

- `grafos/`      leitura/escrita GML e geradores de instâncias (aleatórias e Stacked Book)
- `coloracao/`   coloração gulosa e verificação de coloração de Grundy
- `limites/`     decomposição de vértices, limite superior Stair Factor e limite inferior LB-DEG
- `cliques/`     Bron–Kerbosch (cliques maximais / conjuntos independentes maximais)
- `exatos/`      algoritmos exatos: Programação Dinâmica (PDG) e Branch and Bound (BBG)
- `pli/`         formulações de Programação Linear Inteira (F_N, F_DELTA, F_STF, F_PROD)
                 e modelos auxiliares de coloração própria

Convenção: os vértices em Julia são `1:n`; o vértice `v` corresponde ao rótulo `v-1`
das instâncias originais em Python.
"""
module NumberGrundy

using Graphs
using JuMP
using HiGHS
using Random

# Grafos
include("grafos/io_gml.jl")
include("grafos/geradores.jl")

# Coloração
include("coloracao/gulosa.jl")

# Limites
include("limites/decomposicao.jl")
include("limites/stair_factor.jl")
include("limites/lb_deg.jl")

# Cliques
include("cliques/bron_kerbosch.jl")

# Algoritmos exatos
include("exatos/resultado.jl")
include("exatos/programacao_dinamica.jl")
include("exatos/branch_and_bound.jl")

# Programação Linear Inteira
include("pli/formulacao_grundy.jl")
include("pli/modelos_coloracao.jl")

export
    # grafos
    ler_gml, escrever_gml, grafo_aleatorio, grafo_stacked_book,
    # coloração
    cor_minima, coloracao_gulosa, eh_coloracao_grundy,
    # limites
    grau_maximo, decomposicao_vertices, stair_factor, lb_deg,
    # cliques
    cliques_maximais, conjuntos_independentes_maximais,
    # exatos
    Resultado, pdg, bbg,
    # pli
    pli_grundy, modelo_grundy, limite_formulacao, FORMULACOES,
    pli_numero_cromatico
end
