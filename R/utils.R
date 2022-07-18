#'  Ensure matrix
#'
#' Verify if `data` is a matrix. If not, it is converted into
#' a matrix if possible. Else, an error is thrown.
#' @param data The data
#' @return The data as a matrix
ensure_matrix <- function(data) {
    
    if(!is.matrix(data)){
        tryCatch({
            data = as.matrix(data)
        },
        error = function(){
            message(paste(deparse(substitute(data))," must be a matrix."))
            message(paste(deparse(substitute(data)),"could not be converted into a matrix."))
        })
    }
    
    return(data)
    
}

