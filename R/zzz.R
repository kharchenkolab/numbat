
.onLoad <- function(libname, pkgname){
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
    data.table::setDTthreads(1)
}
