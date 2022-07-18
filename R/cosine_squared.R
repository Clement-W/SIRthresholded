#'  Cosine squared
#'
#' Calculates the quality of the correlation between two vectors
#' (real beta and estimated beta)
#' @param b1 First vector
#' @param b2 Second vector
#' @return The cosine squared between the two vectors
cosine_squared <- function(b1, b2) {
    (matrix(b1, nrow = 1) %*% matrix(b2, ncol = 1)) ^ 2 / ((matrix(b1, nrow = 1)
    %*% matrix(b1, ncol = 1)) * (matrix(b2, nrow = 1) %*% matrix(b2, ncol = 1)))
}