# Thin ccall wrappers around the actual C++ optionspricer implementation
# (julia/csrc/binomial_c_api.cpp), used only by the test suite as the
# ground truth to port against.
module CppReference

const CSRC_DIR = normpath(joinpath(@__DIR__, "..", "csrc"))
const LIB_PATH = joinpath(CSRC_DIR, Sys.isapple() ? "libbinomial_c_api.dylib" : "libbinomial_c_api.so")

function ensure_built()
    isfile(LIB_PATH) && return LIB_PATH
    run(`$(joinpath(CSRC_DIR, "build.sh"))`)
    isfile(LIB_PATH) || error("build.sh ran but $LIB_PATH was not produced")
    return LIB_PATH
end

function binomial_price(S, K, r, q, sigma, T, type::Int, style::Int, steps::Integer)
    ccall((:binomial_price_c, LIB_PATH), Float64,
          (Float64, Float64, Float64, Float64, Float64, Float64, Cint, Cint, Cuint),
          S, K, r, q, sigma, T, type, style, steps)
end

function black_scholes_price(S, K, r, q, sigma, T, type::Int)
    ccall((:black_scholes_price_c, LIB_PATH), Float64,
          (Float64, Float64, Float64, Float64, Float64, Float64, Cint),
          S, K, r, q, sigma, T, type)
end

end # module CppReference
