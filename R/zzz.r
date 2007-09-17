.onLoad <- function (libname, pkgname)
{
  #library.dynam("GlobalAncova", pkgname, libname)
  #dll <- file.path(libname, paste("GlobalAncova", .Platform$dynlib.ext, sep=""))
  #dyn.load(dll)   
  require(methods)
  #require("MASS")
  #require("haplo.stats")
}

.onAttach <- function(lib, pkg) {
    if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
      addVigs2WinMenu("GlobalAncova")
    }
}

.onUnload <- function( libpath ) {
  library.dynam.unload("GlobalAncova", libpath)
}



