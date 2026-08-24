module OptionsPricer

export OptionType, Call, Put
export ExerciseStyle, European, American
export binomial_price

abstract type OptionType end
struct Call <: OptionType end
struct Put <: OptionType end

abstract type ExerciseStyle end
struct European <: ExerciseStyle end
struct American <: ExerciseStyle end

payoff(::Call, S, K) = max(S - K, zero(S))
payoff(::Put, S, K) = max(K - S, zero(S))

# European nodes never compare against immediate exercise, so this method
# never computes the node's underlying price at all.
function relax!(values::AbstractVector, ::European, type::OptionType,
                 S, K, u, d, p, discount, step)
    @inbounds for i in 0:(step - 1)
        values[i + 1] = discount * (p * values[i + 2] + (1 - p) * values[i + 1])
    end
    return values
end

function relax!(values::AbstractVector, ::American, type::OptionType,
                 S, K, u, d, p, discount, step)
    @inbounds for i in 0:(step - 1)
        continuation = discount * (p * values[i + 2] + (1 - p) * values[i + 1])
        S_node = S * u^i * d^(step - 1 - i)
        values[i + 1] = max(continuation, payoff(type, S_node, K))
    end
    return values
end

"""
    binomial_price(S, K, r, q, sigma, T, steps, type, style)

Cox-Ross-Rubinstein binomial tree price of a vanilla option. `r` and `q` are
continuously-compounded (rate, dividend yield); `r` may be negative.
`type` is `Call()`/`Put()`, `style` is `European()`/`American()`.
"""
function binomial_price(S::T, K::T, r::T, q::T, sigma::T, Tmat::T,
                         steps::Integer, type::OptionType,
                         style::ExerciseStyle) where {T<:AbstractFloat}
    steps > 0 || throw(ArgumentError("steps must be positive"))

    dt = Tmat / steps
    u = exp(sigma * sqrt(dt))
    d = inv(u)
    p = (exp((r - q) * dt) - d) / (u - d)
    discount = exp(-r * dt)

    values = Vector{T}(undef, steps + 1)
    @inbounds for i in 0:steps
        S_T = S * u^i * d^(steps - i)
        values[i + 1] = payoff(type, S_T, K)
    end

    for step in steps:-1:1
        relax!(values, style, type, S, K, u, d, p, discount, step)
    end

    return values[1]
end

function binomial_price(S::Real, K::Real, r::Real, q::Real, sigma::Real, Tmat::Real,
                         steps::Integer, type::OptionType, style::ExerciseStyle)
    F = float(promote_type(typeof(S), typeof(K), typeof(r), typeof(q),
                            typeof(sigma), typeof(Tmat)))
    return binomial_price(F(S), F(K), F(r), F(q), F(sigma), F(Tmat), steps, type, style)
end

end # module OptionsPricer
