mcmctreefiles <- grep("mcmctree", dir(), value = T)
for(i in 1:length(mcmctreefiles)){
      setwd(mcmctreefiles[i])
      ctl <- readLines("dating_in.ctl")
      ctl[grep("print", ctl)] <- "print = 1"
      system("rm dating_in.ctl")
      writeLines(ctl, con = "dating_in.ctl")
      setwd("..")
}