# Function to compute the orthonormal matrix from a flat vector using QR decomposition
orthonormalize <- function(A,q) {
  qr.Q(qr(matrix(A, ncol = q)))
}