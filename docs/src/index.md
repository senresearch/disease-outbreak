# DiseaseOutbreak.jl Documentation
```@meta
CurrentModule = DiseaseOutbreak
```
```@contents
```
## Functions
```@docs
<!-- writing function names without parameters will include all the functions that has that name, even the ones from other module. There are two ways to solve this problem. Either include the type in the signature with functionName(::T), or or declare the specific modules that makedocs should include with
makedocs(
    # options
    modules = [MyModule]
)  -->
change(s::Vector{Float64},d::SIRX)

<!-- packageuri(pkgname::String)
tokenisdefined()
token()
githubauth()
repoid(package_name::String)
report(package_name::String, title::String, body::String) -->
```

## Index
```@index
```

```@autodocs
Modules = [DiseaseOutbreak]
```
