'.onAttach' <- function(lib, pkg="bptwin")
  {    
    desc <- packageDescription(pkg)
    packageStartupMessage("\nLoading '", desc$Package, "' package...\n",
                          "Version    : ", desc$Version, "\n",
                          "Overview: help(package=", desc$Package, ")\n");
  }
