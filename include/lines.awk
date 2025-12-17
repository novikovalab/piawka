@namespace "calc"

BEGIN{
  stats::add_stat("lines", "number of lines used in calculation", 0, "")
}

function increment_lines(i){
  thisnum[i]["lines"]=1
  thisden[i]["lines"]=1
}

function finalize_lines(i){
  den[i]["lines"]=1
}
