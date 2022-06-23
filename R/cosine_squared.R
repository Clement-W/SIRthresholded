#'  Cosine squared
#'
#' Calculates the quality of the correlation between two vectors
#' (real beta and estimated beta)
#' @param b1 First vector
#' @param b2 Second vector
#' @return The cosine squared of the two vectors
#' @examples
#' b1 = c(-0.22,-0.31,-0.27)
#' b2 = c(1,1,1)
#' cosine_squared(b1,b2)
cosine_squared <- function(b1, b2) {
    (matrix(b1, nrow = 1) %*% matrix(b2, ncol = 1)) ^ 2 / ((matrix(b1, nrow = 1)
    %*% matrix(b1, ncol = 1)) * (matrix(b2, nrow = 1) %*% matrix(b2, ncol = 1)))
}