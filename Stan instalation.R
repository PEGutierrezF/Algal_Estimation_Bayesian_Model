
#https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows
# https://discourse.mc-stan.org/t/rstan-wont-compile-with-r3-6/9105/4

install.packages("rtools40")
pkgbuild::has_build_tools(debug = TRUE)

remove.packages("rstan")
if (file.exists(".RData")) file.remove(".RData")

file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
M <- file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
file.edit(M)

Sys.setenv(MAKEFLAGS = paste0("-j", parallel::detectCores()))
install.packages("pkgbuild", INSTALL_opts = "--no-multiarch")
cat("Rtools version 4.0.0", file = file.path("C:", "rtools40", "VERSION.txt"), sep = "\n")
cat("Rtools 4.0", file = file.path("C:", "rtools40", "Rtools.txt"), sep = "\n")
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE, INSTALL_opts = "--no-multiarch")
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)



Sys.setenv(MAKEFLAGS = "-j4") # four cores used
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)


library(usethis)
library(devtools)


