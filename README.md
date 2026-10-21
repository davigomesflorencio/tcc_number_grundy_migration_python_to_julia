# NumberGrundy.jl — Algoritmos para o Problema do Número de Grundy

Migração para **Julia** do repositório do TCC
[davigomesflorencio/number_grundy](https://github.com/davigomesflorencio/number_grundy)
(Ciência da Computação, UFC Campus Quixadá), originalmente em Python
(NetworkX + Google OR-Tools).

O número de Grundy (ou número cromático guloso) de um grafo é a maior quantidade de
cores que o Algoritmo Guloso de Coloração consegue usar, considerando todas as ordens
dos vértices. Determinar esse número é NP-completo.

## Correspondência Python → Julia

| Algoritmo (TCC)                    | Python                              | Julia                                  | Arquivo |
|-----------------------------------|-------------------------------------|----------------------------------------|---------|
| Coloração gulosa                   | `coloracaoGulosa`, `corMinima`      | `coloracao_gulosa`, `cor_minima`       | `src/coloracao/gulosa.jl` |
| Verificação de coloração de Grundy | `isGrundyColoring`                  | `eh_coloracao_grundy`                  | `src/coloracao/gulosa.jl` |
| Decomposição de vértices           | `ordemCrescenteGrau`, `vertex_decomposition` | `decomposicao_vertices`       | `src/limites/decomposicao.jl` |
| Limite superior **Stair Factor**   | `stairfactor`                       | `stair_factor`                         | `src/limites/stair_factor.jl` |
| Limite inferior **LB-DEG**         | `sgNumberGrundy`                    | `lb_deg`                               | `src/limites/lb_deg.jl` |
| Bron–Kerbosch                      | `findCliques`, `findCliques3`       | `cliques_maximais`, `conjuntos_independentes_maximais` | `src/cliques/bron_kerbosch.jl` |
| **PDG** — Programação Dinâmica     | `lawler`                            | `pdg`                                  | `src/exatos/programacao_dinamica.jl` |
| **BBG** — Branch and Bound         | `BAB`, `BABNumberGrundy`            | `bbg`                                  | `src/exatos/branch_and_bound.jl` |
| **PLI** F_N / F_DELTA / F_STF / F_PROD | `model_grundy(G, ub)`           | `pli_grundy(G, :n / :delta / :stf / :prod)` | `src/pli/formulacao_grundy.jl` |
| PLI de coloração própria (χ)       | modelos do notebook `Original/`     | `pli_numero_cromatico`                 | `src/pli/modelos_coloracao.jl` |
| Grafos aleatórios G(n,p)           | `generateGraph`                     | `grafo_aleatorio`                      | `src/grafos/geradores.jl` |
| Grafos Stacked Book                | `generateStackedBookGraph`          | `grafo_stacked_book`                   | `src/grafos/geradores.jl` |
| Leitura/escrita GML                | `nx.read_gml` / `nx.write_gml`      | `ler_gml` / `escrever_gml`             | `src/grafos/io_gml.jl` |

### Formulações PLI

Todas usam o mesmo modelo (variáveis `x[v,c]`, `z[c]` binárias) e diferem no número
de cores `|C|`:

| Formulação | `|C|`            | Arquivo de resultados do TCC |
|------------|------------------|------------------------------|
| `:n`       | n                | `dados_random/PLI_N.csv` |
| `:delta`   | Δ(G) + 1         | `PLI_DELTA.csv` |
| `:stf`     | Stair Factor     | `PLI_STF.csv` |
| `:prod`    | 6                | `dados/PLI_6.csv` (Stacked Book) |

O solver padrão é o **HiGHS** (open source, instalado automaticamente). O TCC usou o
**SCIP**; para usá-lo, adicione `SCIP` ao ambiente e passe o otimizador:

```julia
using SCIP
pli_grundy(g, :stf; otimizador = SCIP.Optimizer)
```

## Estrutura

```
src/
  NumberGrundy.jl              módulo principal
  grafos/                      io_gml.jl, geradores.jl
  coloracao/                   gulosa.jl
  limites/                     decomposicao.jl, stair_factor.jl, lb_deg.jl
  cliques/                     bron_kerbosch.jl
  exatos/                      resultado.jl, programacao_dinamica.jl, branch_and_bound.jl
  pli/                         formulacao_grundy.jl, modelos_coloracao.jl
experimentos/                  scripts dos testes computacionais
test/runtests.jl               testes (inclui comparação com os resultados do TCC)
instancias/                    instâncias originais (random/*.gml, classe/*)
dados/, dados_random/, tempo/, TabelasTCC/   resultados originais do TCC (Python)
python_original/               notebooks Python originais, para referência
```

## Instalação

Requer Julia ≥ 1.10.

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

## Uso

```julia
using NumberGrundy, Graphs

g = ler_gml("instancias/classe/classe_graph_3_4")   # ou grafo_stacked_book(3, 4)

stair_factor(g)                 # limite superior
lb_deg(g)                       # limite inferior (valor, coloração)

pdg(g)                          # Programação Dinâmica (n ≤ 30)
bbg(g, stair_factor(g))         # Branch and Bound
pli_grundy(g, :stf)             # PLI F_STF
pli_grundy(g, 6; tempo_limite = 60)
```

Todos os algoritmos exatos retornam um `Resultado` com os campos `grundy`, `otimo`,
`tempo` e `coloracao`. Com `tempo_limite` (padrão 900 s, como no TCC) o algoritmo
devolve o melhor valor encontrado e `otimo = false` — substitui o `signal.alarm` e a
variável global `config.k` do Python.

## Experimentos

```bash
# Grafos aleatórios (n ∈ {10,15,20,25}, p ∈ {0.1,...,0.9}, 5 repetições)
julia --project=. experimentos/experimento_random.jl PLI_DELTA [tempo_limite] [saida.csv]

# Stacked Book (n ∈ 2..9, m ∈ 3..9)
julia --project=. experimentos/experimento_stacked_book.jl BB_STF [tempo_limite] [saida.csv]

# Gerar novas instâncias (em instancias_geradas/)
julia --project=. experimentos/gerar_instancias.jl
```

Algoritmos: `PDG`, `BB_DELTA`, `BB_STF`, `BB_6`, `PLI_N`, `PLI_DELTA`, `PLI_STF`, `PLI_6`.
Os CSVs são gravados em `resultados/` com as mesmas colunas dos arquivos do TCC.

## Testes

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

Os testes conferem que, nas instâncias do TCC, o gerador Stacked Book, Δ+1, Stair
Factor e o número de Grundy (PDG, BBG e PLI) reproduzem os valores dos CSVs originais.

## Diferenças em relação ao Python

- Vértices são `1:n` (o rótulo `k` das instâncias vira o vértice `k+1`); "sem cor" é `0`.
- PDG usa conjuntos em bitmask em vez de strings como chave de dicionário, e o
  Bron–Kerbosch em bitmask com pivô de Tomita (mesmos conjuntos, muito mais rápido).
- BBG colore no próprio vetor com backtracking em vez de copiar listas a cada nó;
  a ordem de ramificação, as cores testadas e a poda são as mesmas.
- `lb_deg`: o `sgNumberGrundy` original tinha um erro de índice (`c[index]` no lugar
  de `c[y]`) e podia gerar colorações que não são de Grundy; a versão Julia aplica a
  coloração gulosa em ordem decrescente de grau, o que sempre gera um limite válido.
- A decomposição de vértices foi traduzida literalmente para reproduzir o Stair
  Factor do TCC (o comentário em `decomposicao.jl` explica o que ela calcula).
