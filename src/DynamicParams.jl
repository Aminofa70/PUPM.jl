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