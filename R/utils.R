# Verify if `data` is a matrix. If not, it is converted into
# a matrix if possible. Else, an error is thrown.
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

# Calculates the quality of the correlation between two vectors
cosine_squared <- function(b1, b2) {
    (matrix(b1, nrow = 1) %*% matrix(b2, ncol = 1)) ^ 2 / ((matrix(b1, nrow = 1)
    %*% matrix(b1, ncol = 1)) * (matrix(b2, nrow = 1) %*% matrix(b2, ncol = 1)))
}

