# Limite inferior LB-DEG: heurística gulosa guiada pelo grau.

"""
    lb_deg(g) -> (lb, coloracao)

Limite inferior para o número de Grundy. Os vértices são ordenados por grau
decrescente e a coloração gulosa é aplicada nessa ordem; toda coloração gulosa é
uma coloração de Grundy, portanto o número de cores usado é um limite inferior
válido.

Observação: no notebook original (`sgNumberGrundy`) a heurística recoloria os
vizinhos de cada vértice e podia gerar colorações que não são de Grundy (havia
um erro de índice `c[index]` no lugar de `c[y]`). Aqui usa-se a versão correta.
"""
function lb_deg(g::SimpleGraph)
    ordem = sort(collect(vertices(g)); by=v -> degree(g, v), rev=true)
    return coloracao_gulosa(g, ordem)
end
