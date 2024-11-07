
Normalize = function(B) {
  # B (p x q)-matrix
  B<-matrix(B,ncol=q)
  if (q==1) {
    sg<-B[1,1]/abs(B[1,1])
    return(sg*B/norm(B,"2"))
  }
  Brrefortho=gramSchmidt(t(rref(t(B))))$Q
  
  return(Brrefortho)
}
