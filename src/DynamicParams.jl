# src/DynamicParams.jl
mutable struct DynamicParams
    parameters::Dict{Symbol, Any}
    function DynamicParams()
        new(Dict{Symbol, Any}())
    end
end

# Overriding getproperty
function Base.getproperty(obj::DynamicParams, name::Symbol)
    if name == :parameters || haskey(obj.parameters, name)
        return name == :parameters ? getfield(obj, :parameters) : obj.parameters[name]
    else
        throw(ArgumentError("Property '$name' does not exist."))
    end
end

# Overriding setproperty!
function Base.setproperty!(obj::DynamicParams, name::Symbol, value)
    if name == :parameters
        throw(ArgumentError("Cannot modify the `parameters` field directly."))
    else
        obj.parameters[name] = value
    end
end

function my_mean(data::Vector{<:Number})
    total = 0.0
    for x in data
        total += x
    end
    return total / length(data)
end

function my_std(data::Vector{<:Number})
    μ = my_mean(data)
    sumsq = 0.0
    for x in data
        sumsq += (x - μ)^2
    end
    return sqrt(sumsq / (length(data) - 1))  # sample standard deviation
end