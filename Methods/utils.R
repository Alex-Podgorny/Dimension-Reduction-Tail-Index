

# Function to compute the orthonormal matrix from a flat vector using QR decomposition
orthonormalize <- function(A,q) {
  qr.Q(qr(matrix(A, ncol = q)))
}

# Function to compute the canonical basis
Normalize = function(B,q) {
  # B (p x q)-matrix
  B<-matrix(B,ncol=q)
  if (q==1) {
    sg<-B[1,1]/abs(B[1,1])
    return(sg*B/norm(B,"2"))
  }
  Brrefortho=orthonormalize(t(rref(t(B))),q)
  
  return(Brrefortho)
}