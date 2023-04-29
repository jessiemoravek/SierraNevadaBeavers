# evenspace creates points along lines with a set distance
evenspace <- function(xy, sep, start=0){
  
  dx <- c(0,diff(xy[,1]))
  dy <- c(0,diff(xy[,2]))
  dseg <- sqrt(dx^2+dy^2)
  dtotal <- cumsum(dseg)
  
  linelength = sum(dseg)
  
  pos = seq(start,linelength, by=sep)
  pos = pos[1:(length(pos)-1)]  # to get rid of the last point because it is enclosed
  
  whichseg = unlist(lapply(pos, function(x){sum(dtotal<=x)}))  # return index of elements that the accumulative length < segmentation length
  #whichseg = whichseg[1:(length(whichseg)-1)]
  
  pos = data.frame(pos=pos,whichseg=whichseg,
                   x0=xy[whichseg,1],
                   y0=xy[whichseg,2],
                   dseg = dseg[whichseg+1],
                   dtotal = dtotal[whichseg],
                   x1=xy[whichseg+1,1],
                   y1=xy[whichseg+1,2]
  )
  
  pos$further =  pos$pos - pos$dtotal  #the segmentation might not be exact the integers of set distance.
  pos$f = pos$further/pos$dseg
  pos$x = pos$x0 + pos$f * (pos$x1-pos$x0) # to make sure the segmentation x coord is exactly as required
  pos$y = pos$y0 + pos$f * (pos$y1-pos$y0)# to make sure the segmentation y coord is exactly as required
  
  pos$theta = atan2(pos$y0-pos$y1,pos$x0-pos$x1)
  
  return(pos[,c("x","y","x0","y0","x1","y1","theta")])
}

