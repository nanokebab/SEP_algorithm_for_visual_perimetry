makeStimHelper <-
  function(db,n, x, y) {
    # returns a function of (db,n)
    ff <- function(db, n)
      db + n
    body(ff) <- substitute({
      s <- list(
        x = x, y = y, level = dbTocd(db), size = 0.43, color = "white",
        duration = 200, responseWindow = 1500
      )
      class(s) <- "opiStaticStimulus"
      return(s)
    }
    , list(x = x,y = y))
    return(ff)
  }