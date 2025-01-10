"""
function to remove vtu file 
```
remove_vtk_files(directory::String)
```
"""
function remove_vtk_files(directory::String)
    files = readdir(directory)
    for file in files
        if endswith(file, ".vtu")
            filepath = joinpath(directory, file)
            rm(filepath, force = true)
        end
    end
end


