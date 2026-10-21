# Tipo de retorno comum dos algoritmos exatos e das formulações PLI.

"""
    Resultado

- `grundy`: melhor valor encontrado (ótimo se `otimo == true`)
- `otimo`: `true` se o algoritmo terminou provando a otimalidade
- `tempo`: tempo de execução em segundos
- `coloracao`: coloração de Grundy que atinge `grundy` (vazia se indisponível)

No Python o limite de tempo era feito com `signal.alarm(900)` e o melhor valor
parcial era guardado na variável global `config.k`; aqui cada algoritmo recebe
`tempo_limite` e devolve o melhor valor parcial com `otimo = false`.
"""
struct Resultado
    grundy::Int
    otimo::Bool
    tempo::Float64
    coloracao::Vector{Int}
end

Base.show(io::IO, r::Resultado) =
    print(io, "Resultado(grundy=$(r.grundy), otimo=$(r.otimo), tempo=$(round(r.tempo; digits=4))s)")

struct TempoEsgotado <: Exception end

@inline _verificar_prazo(prazo::Float64) = time() > prazo && throw(TempoEsgotado())
