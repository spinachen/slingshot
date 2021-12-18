sling <- function(x,y,p,l){
  reticulate::source_python("./slingshot/R/sling.py")
  x = scale(x)
  result = Sling(x,y,p,l)
  return(list(w=result[1],iter=result[2],obj=result[3]))
}