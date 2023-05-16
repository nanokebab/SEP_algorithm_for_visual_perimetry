dynamic.stop <- function (state) 
{
  return(state$finished != "Not")
}
environment(dynamic.stop) <- asNamespace('OPI')