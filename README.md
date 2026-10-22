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
  v2/                          branch_and_bound.jl (BBG v2, ver seção "v2")
  v3/                          k_atomos.jl, limite_vertices.jl (BBG v3, ver seção "v3")
  v4/                          limite_arestas.jl, k_atomos_arestas.jl (BBG v4, ver seção "v4")
experimentos/                  scripts dos testes computacionais
  v2/                          comparar_bbg.jl (BBG original × BBG v2)
  v3/                          avaliar_bbg_v3.jl (BBG v3 × melhores valores do TCC)
  v4/                          avaliar_bbg_v4.jl (BBG v4 × v3 × valores do TCC)
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

## v2 — Branch and Bound melhorado (`src/v2/`)

`bbg_v2(g, ub)` mantém a árvore de busca do BBG do TCC e acrescenta:

1. **Poda válida.** O BBG original podava quando a maior cor usada *até o momento*
   era menor que o LB, mas os vértices restantes ainda podem receber cores maiores.
   Por isso ele **não é exato**: em 2000 grafos aleatórios pequenos devolveu um valor
   menor que o ótimo em 33. A v2 poda quando
   `max(cor_maxima, max_{u não colorido} min(d(u)+1, UB)) ≤ LB`.
2. **Parada antecipada** quando LB = UB.
3. **LB inicial heurístico** (LB-DEG e coloração gulosa na ordem de decomposição) em
   vez de 0.
4. **Verificação incremental da propriedade de Grundy.** Um vértice de cor `c` precisa
   de vizinhos com todas as cores `1..c-1`, e cada vizinho sem cor fornece no máximo
   uma. O ramo é podado assim que `faltam(w) > vizinhos_sem_cor(w)`, a um custo de
   O(Δ) por atribuição. Toda folha alcançada já é uma coloração de Grundy.

```julia
bbg_v2(g, grau_maximo(g) + 1)      # exato
```

```bash
# BBG original × BBG v2 nas instâncias do TCC (UB = DELTA, STF ou 6; limite em segundos)
julia --project=. experimentos/v2/comparar_bbg.jl STF 30
```

Algoritmos `BB2_DELTA`, `BB2_STF` e `BB2_6` também estão disponíveis nos scripts de
experimentos.

> **Atenção — Stair Factor:** a função `stair_factor` reproduz fielmente o cálculo do
> TCC, mas **não é um limite superior válido em todo grafo**. Exemplo: arestas
> `(1,3),(1,4),(2,4),(2,5),(3,5),(3,6),(5,6)` têm Γ = 4 e `stair_factor` = 3. Nas
> instâncias do TCC isso não ocorre. Para garantir o ótimo use UB = Δ+1 ou
> `limite_superior_vertices(g)` (v3).

## v3 — Busca de testemunhas / k-átomos (`src/v3/`)

`bbg_v3(g)` não colore o grafo inteiro. Ele usa o fato de que

> Γ(G) ≥ k ⇔ existe uma **coloração parcial de Grundy** (coloração própria de um
> subconjunto S em que todo vértice de cor c tem, em S, vizinhos com as cores
> 1..c-1) com algum vértice de cor k.

Essa coloração parcial sempre se estende a uma coloração de Grundy completa pelo
algoritmo guloso (vértices de S primeiro, em ordem crescente de cor). A busca
escolhe uma raiz com cor k e, a cada passo, a demanda mais restrita (vértice `w` de
cor `c` sem vizinho de cor `j < c`), ramificando sobre qual vizinho de `w` recebe a
cor `j`. A testemunha tem no máximo **2^(k-1)** vértices (o tamanho de um k-átomo),
qualquer que seja n.

- **Teste "Γ ≥ k?"**: O(n · k · W · Δ^W), com W = min(n, 2^(k-1)) — n raízes, árvore
  de profundidade ≤ W com ramificação ≤ Δ, e O(W · k · Δ) por nó (percorre os
  vértices coloridos, suas cores faltantes e seus vizinhos). Polinomial para k fixo:
  O(n · Δ^(2^(k-1))) ⊆ O(n^(2^(k-1)+1)).
- **Algoritmo**: testa k = LB+1, LB+2, … até a primeira resposta negativa. Cada
  resposta positiva é completada gulosamente e pode saltar vários valores de k.
- **Total**: dominado pelo teste negativo k = Γ+1:
  **O(n · Γ · W · Δ^W)** com W = min(n, 2^Γ). Memória O(k·n + m + Δ·W).

| Versão | Pior caso | Expoente depende de |
|--------|-----------|---------------------|
| v1 | O((n+m)·(Δ+1)ⁿ) | n |
| v2 | O(Δ·(Δ+1)ⁿ) | n |
| v3 | O(n·Γ·W·Δ^W), W = min(n, 2^Γ) | Γ (limitado por n) |

O v3 também traz um **limite superior por vértice**, `limites_por_vertice(g)`:
partindo de `b(v) = d(v)+1`, itera-se até o ponto fixo
`b(v) = 1 + max{t : v tem vizinhos u_1..u_t com b(u_i) ≥ i}`. O valor
`limite_superior_vertices(g) = max b(v)` é **sempre** um limite superior válido para
Γ(G), ao contrário do `stair_factor` do TCC. Também restringe quais vértices podem
receber cada cor durante a busca.

```julia
bbg_v3(g)                          # exato, UB = limite_superior_vertices(g)
bbg_v3(g; tempo_limite = 60)
```

```bash
# v3 nas instâncias do TCC, comparado com o melhor valor obtido no TCC
julia --project=. experimentos/v3/avaliar_bbg_v3.jl 30
```

## v4 — Limite superior por aresta (`src/v4/`)

`bbg_v4(g)` é a busca de k-átomos da v3 com um limite mais forte. Para cada aresta
orientada u→v, a mensagem m(u→v) limita a cor de `u` **quando `v` tem cor maior**
(isto é, quando `u` é testemunha de `v`). Nesse caso `v` não pode ser testemunha de
`u`, então só os outros vizinhos contam:

```
m(u→v) = 1 + T({ m(x→u) : x ∈ N(u) \ {v} })        (ponto fixo, partindo de d(u))
b(v)   = 1 + T({ m(u→v) : u ∈ N(v) })
T(β)   = max{ t : β em ordem decrescente, β_i ≥ t − i + 1 para i = 1..t }
```

`limite_superior_arestas(g) = max b(v)` é sempre válido e nunca maior que o
`limite_superior_vertices` da v3. Na busca, um vizinho `u` só é candidato a
testemunha de cor `j` para `w` se m(u→w) ≥ j.

| Grafo | Δ+1 | Stair Factor | v3 (vértice) | **v4 (aresta)** | Γ |
|-------|-----|--------------|--------------|-----------------|---|
| Estrela K₁,₆ | 7 | — | 3 | **2** | 2 |
| Stacked Book 2×9 | 11 | 4 | 5 | **4** | 4 |
| Stacked Book 9×9 | 12 | 12 | 7 | **6** | 6 |

- Em **todas as 56 instâncias Stacked Book** do TCC, o limite por aresta é igual a
  Γ. Nessa classe a v4 só precisa encontrar a testemunha, sem refutar Γ+1.
- Em 1650 grafos aleatórios pequenos, a distância média até Γ caiu de 0,99 (v3) para
  0,59 (v4).
- Para Stacked Book, o limite 6 do TCC (F_PROD) tem uma prova direta: folhas têm
  grau ≤ 3; se o centro tem cor ≥ 4, as folhas dele só usam as duas folhas vizinhas
  como testemunhas (cor ≤ 3); então um centro tem cor ≤ 1 + 3 + 2 = 6. O limite por
  aresta reproduz esse argumento automaticamente.

```julia
limite_superior_arestas(g)        # UB válido
bbg_v4(g; tempo_limite = 60)      # exato
pli_grundy(g, limite_superior_arestas(g))   # PLI com o UB por aresta
```

```bash
julia --project=. experimentos/v4/avaliar_bbg_v4.jl 30
```

Ainda não implementado na v4: tabela de transposição, poda de Hall, quebra de
simetria nas raízes e LB por busca local.

## Resultados: BBG v1 × v2 × v3 × v4

Instâncias resolvidas (otimalidade provada) nas 156 instâncias do TCC.

**Configuração:** 30 s por instância e por algoritmo; Julia 1.13, HiGHS, Windows 11,
uma thread. v1 e v2 usam UB = Stair Factor (como o BB_STF do TCC); v3 e v4 usam os
próprios limites. A comparação v1 × v2 rodou ao mesmo tempo que as avaliações da v3
e da v4, então os tempos medidos podem estar um pouco piores que o normal. A linha
"TCC" conta uma instância como resolvida se **algum** algoritmo do TCC (PD, BB ou
PLI, com SCIP) terminou em até 900 s.

| Algoritmo | Stacked Book (56) | Aleatórias (100) | **Total (156)** |
|-----------|-------------------|------------------|-----------------|
| v1 — BBG do TCC | 11 | 42 | **53** |
| v2 — poda válida, LB heurístico, verificação incremental | 28 | 55 | **83** |
| v3 — k-átomos, limite por vértice | 51 | 67 | **118** |
| **v4 — k-átomos, limite por aresta** | **56** | **69** | **125** |
| TCC (melhor entre PD, BB e PLI, 900 s) | 56 | 72 | 128 |

- De v1 para v4, as instâncias resolvidas passam de 53 para 125, com o mesmo limite de
  30 s. A v4 fica a 3 instâncias do TCC usando 1/30 do tempo.
- A v4 resolve **todas** as Stacked Book, cada uma em milissegundos: o limite por
  aresta é igual a Γ, então basta encontrar a testemunha.

  | Stacked Book | v3 | v4 |
  |--------------|----|----|
  | 5×9 | 5, estourou 30 s | 5, ótimo, 0,0 s |
  | 6×9 a 9×9 | 5, estourou 30 s | 6, ótimo, 0,0 s |

- Em 80 instâncias, a v4 encontra Γ igual ao seu limite superior, e a prova de ótimo
  sai sem refutar Γ+1.
- Em 5 aleatórias não resolvidas, a v4 encontra valores maiores que a v3 (por
  exemplo, `random_4_25_0.7`: 13 → 15).
- Nas instâncias resolvidas, v1 e v2 deram o mesmo valor que a v4. Os problemas do v1
  (poda que pode descartar o ótimo) e do Stair Factor (que pode ficar abaixo de Γ)
  aparecem em grafos aleatórios pequenos, mas não nas instâncias do TCC.
- **Onde a v4 ainda não resolve:** 31 aleatórias densas (p = 0,9; p = 0,7 com n ≥ 20;
  p = 0,5 com n = 25, além de `random_0_20_0.5`). Nelas o limite por aresta fica perto
  de n (por exemplo, 24 ou 25 com Γ ≈ 16), e a busca trava na refutação de Γ+1.
  Possíveis próximos passos: tabela de transposição, poda de Hall e usar a PLI quando
  a busca empacar.

Para reproduzir (os CSVs vão para `resultados/`, que fica fora do Git):

```bash
julia --project=. experimentos/v2/comparar_bbg.jl STF 30     # v1 × v2
julia --project=. experimentos/v3/avaliar_bbg_v3.jl 30       # v3
julia --project=. experimentos/v4/avaliar_bbg_v4.jl 30       # v4 (lê o CSV da v3)
```

## Atenção: valores do TCC em execuções que estouraram 900 s

Nos CSVs originais (`dados/`, `dados_random/`), as linhas com `tempo = 900` **não
são confiáveis**. No laço de experimentos do Python, ao estourar o tempo
(`signal.alarm(900)`) era gravado o valor de `config.k`, uma variável global que:

- só era zerada **depois de uma execução bem-sucedida** (nos grafos aleatórios), e
  por isso podia carregar o valor de uma instância anterior;
- não era atualizada pelo `model_grundy` (PLI).

Evidências (todas em linhas com `tempo = 900`), confirmadas por PD, PLI com n cores,
BBG v2 e BBG v3 em Julia, e pelos PDG/PLI do próprio TCC:

| Instância | Valor no CSV | Γ correto |
|-----------|--------------|-----------|
| Stacked Book 2×8, 2×9, 3×5…3×9 | BB_6 = 5, BB_DELTA = 5 (2×8, 2×9) | 4 |
| Aleatórias n = 15 (p ≥ 0,3) | BB_DELTA = 9 em 17 grafos seguidos | 6, 7 ou 8 |
| Aleatórias n = 25, p = 0,1 | PDG = 5 | 4 |

Comparações de qualidade feitas a partir dessas linhas, como "o BBG supera o PBG
na qualidade das soluções", devem ser revistas. O script
`experimentos/v3/avaliar_bbg_v3.jl` só usa os valores de execuções concluídas e
mostra os demais na coluna informativa `grundy_tcc_estourado`.

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
