"""
Speed of light is set to c = 1 (natural units)
"""

const c = 1.0

"""
Apply Lorentz transformation to coordinates (x,t) with velocity parameter β
Returns (x', t') in the new frame
"""
function lorentz(x::Real, t::Real, β::Real)::Tuple{Float64, Float64}
    abs(β) < 1.0 || throw(ArgumentError("Absolute β must be less than 1"))
    γ = inv(sqrt(1.0 - β^2)) 
    x_prime = γ * (x -  β * c * t)
    t_prime = γ * (t - (β/c) * x)
    
    return(x_prime, t_prime)
end

function inverse_lorentz(x_prime::Real, t_prime::Real, β::Real)::Tuple{Float64, Float64}
   abs(β) < 1.0 || throw(ArgumentError("Absolute β must be less than 1"))
   return lorentz(x_prime, t_prime, -β)
end

function test_inverse_lorentz(x::Real, t::Real, β::Real)::Tuple{Float64, Float64, Real}

# Or use the unified function
   (x_prime, t_prime) = lorentz(x, t, β)
   (x_back, t_back) = inverse_lorentz(x_prime, t_prime,  β)
   println((round(x_back;digits=6),round(t_back;digits=6),β))
   return (x_back, t_back, β)

end

"""
    spacetime_invariant(x::Real, t::Real) -> Float64

Calculate the spacetime interval invariant √|c²t² - x²|
"""
function spacetime_invariant(x::Real, t::Real)::Float64
    interval = (c*t)^2 - x^2
    sqrt(abs(interval))
end

"""
Check if the spacetime interval is timelike (|x| < c|t|)
"""
is_timelike(x::Real, t::Real)::Bool = ( x^2 < (c*t)^2 )


"""
Calculate the Lorentz factor γ = 1/√(1-β²)
"""

function show_lorentz(x::Real, t::Real, β::Real)
    
    # Calculate transformation
    (x_prime, t_prime) = lorentz(x, t, β)
    
    # Display results
    println("\n=== Lorentz Transformation ===")
    println("Original coordinates:")
    println("  x = $x")
    println("  t = $t")
    
    println("\nTransformed coordinates:")
    println("  x' = $(round(x_prime; digits=6))")
    println("  t' = $(round(t_prime; digits=6))")
    
    # Calculate and display invariants
    invariant_original = spacetime_invariant(x, t)
    invariant_transformed = spacetime_invariant(x_prime, t_prime)
    
    interval_type = is_timelike(x, t) ? "timelike" : "spacelike"
    quantity_name = is_timelike(x, t) ? "proper time" : "proper distance"
    
    println("\nInvariant quantities:")
    println("  Interval type: $interval_type")
    println("  $quantity_name (original): $(round(invariant_original; digits=6))")
    println("  $quantity_name (transformed): $(round(invariant_transformed; digits=6))")
    
    # Display parameters
    α = sqrt(1.0 - β^2)
    γ = inv(α)
    println("\nTransformation parameters:")
    println("  c = $c")
    println("  β = $β")
    println("  γ = $(round(γ; digits=6))")
    println("  α = 1/γ = $(round(α; digits=6))")
    
    return nothing
end

