.onLoad <- function (libname, pkgname){
  #require(methods)
}

.onAttach <- function(lib, pkg) {
    if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI == "Rgui"){
      addVigs2WinMenu("GlobalAncova")
    }
}

.onUnload <- function( libpath ) {
  library.dynam.unload("GlobalAncova", libpath)
}



