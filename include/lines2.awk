@namespace "calc"

BEGIN{
  stats::add_stat("lines2", "number of lines used in between-pop calculation", "between", "")
}

function increment_lines2(i,j){
  thisnum[i,j]["lines2"]=1
  thisden[i,j]["lines2"]=1
}

function finalize_lines2(i,j){
  den[i,j]["lines2"]=1
}
