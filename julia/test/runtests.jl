using Test

include(joinpath(@__DIR__, "..", "src", "OptionsPricer.jl"))
using .OptionsPricer

include(joinpath(@__DIR__, "cpp_reference.jl"))
using .CppReference

CppReference.ensure_built()

# Mirrors app/main.cpp's build_grid(): moneyness x vol x maturity x rate x
# {call, put}, K=100, q=0. type_code/style_code follow the C API's encoding
# (0 = Call/European, 1 = Put/American) used by binomial_c_api.cpp.
struct Scenario
    S::Float64
    K::Float64
    r::Float64
    q::Float64
    sigma::Float64
    T::Float64
    type::OptionType
    type_code::Int
end

function build_grid()
    K_STRIKE = 100.0
    moneyness = (0.8, 1.0, 1.2)
    vols = (0.20, 0.40)
    maturities = (0.25, 1.0)
    rates = (0.03, -0.02)
    types = ((Call(), 0), (Put(), 1))
    q = 0.0

    grid = Scenario[]
    for m in moneyness, sigma in vols, T in maturities, r in rates, (type, code) in types
        push!(grid, Scenario(m * K_STRIKE, K_STRIKE, r, q, sigma, T, type, code))
    end
    return grid
end

const GRID = build_grid()
const BINOMIAL_STEPS = 1000

@testset "OptionsPricer.jl" begin

    @testset "agrees with C++ binomial_price across 48-scenario grid" begin
        @test length(GRID) == 48
        for sc in GRID
            jl_price = binomial_price(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T,
                                       BINOMIAL_STEPS, sc.type, European())
            cpp_price = CppReference.binomial_price(sc.S, sc.K, sc.r, sc.q, sc.sigma,
                                                      sc.T, sc.type_code, 0, BINOMIAL_STEPS)
            @test isapprox(jl_price, cpp_price; atol=1e-8)
        end
    end

    @testset "converges to Black-Scholes as steps increase" begin
        S, K, r, q, sigma, T = 100.0, 100.0, 0.05, 0.0, 0.20, 1.0
        bs = CppReference.black_scholes_price(S, K, r, q, sigma, T, 0)

        coarse = binomial_price(S, K, r, q, sigma, T, 25, Call(), European())
        fine = binomial_price(S, K, r, q, sigma, T, 5000, Call(), European())

        @test abs(fine - bs) < abs(coarse - bs)
        @test isapprox(fine, bs; atol=1e-2)

        # O(1/N) signature: N * error should stay roughly bounded rather than
        # growing, across a widening range of step counts.
        step_counts = (10, 50, 100, 500, 1000, 5000)
        n_err = [n * abs(binomial_price(S, K, r, q, sigma, T, n, Call(), European()) - bs)
                 for n in step_counts]
        @test maximum(n_err) < 10 * minimum(n_err)
    end

    @testset "American/European exercise consistency" begin
        steps = 500
        for sc in GRID
            euro = binomial_price(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T,
                                   steps, sc.type, European())
            amer = binomial_price(sc.S, sc.K, sc.r, sc.q, sc.sigma, sc.T,
                                   steps, sc.type, American())

            if sc.type isa Put
                @test amer >= euro - 1e-10
            elseif sc.r >= 0.0  # q=0 grid: American call == European call iff r>=0
                @test isapprox(amer, euro; atol=1e-8)
            else
                @test amer >= euro - 1e-10
            end
        end
    end

end
