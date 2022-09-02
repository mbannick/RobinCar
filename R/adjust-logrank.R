adjust.LogRank <- function(model, data){

  result <- ...

  return(
    structure(
      class="LogRankResult",
      list(result=result, settings=model, data=data)
    )
  )
}
