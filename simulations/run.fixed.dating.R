mcmctreefiles <- grep("mcmctree", dir(), value = T)
for(i in 1:length(mcmctreefiles)){
      setwd(mcmctreefiles[i])
      system("mcmctree dating_in.ctl")
      setwd("..")
}