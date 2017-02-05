table <- read.table("10-50kb_segments_molID.txt", header = TRUE, sep='\t')
segment <- seq(1944)
table$segment <- segment

red_table <- read.table("redlabel_positions_molID.txt", header=TRUE, sep='\t')


png("segments_10-50_redlabels.png", height=3000, width=6000, res=300)
plot(c(table[1, "start"], table[1,"stop"]), rep(table[1,"segment"], 2), type='l', col="blue", xlim=c(0,7500000), ylim=c(0,1944))
red_row <- which(red_table$molID == table[1,"molID"])
if (length(red_row) > 0 ) {
  for (row in red_row) {
    points(red_table[row, "bp"], table[i,"segment"], col="red")
  }
}


for (i in c(2:nrow(table))) {
  red_row <- which(red_table$molID == table[i,"molID"])
  if (length(red_row) > 0 ) {
    for (row in red_row) {
      points(red_table[row, "bp"], table[i,"segment"], col="red", pch=20, cex=0.5)
    }
  }
  lines(c(table[i, "start"], table[i,"stop"]), rep(table[i,"segment"], 2), type='l', col="blue")
 
}
abline(v=3900000, col="red")
dev.off()

#############################################################################################

png("segments_10-50_redlabels_zoom.png", height=5000, width=4000, res=300)
plot(c(table[1, "start"], table[1,"stop"]), rep(table[1,"segment"], 2), type='l', col="blue",  xlim=c(3600000,4200000), ylim=c(0,1944))
red_row <- which(red_table$molID == table[1,"molID"])
if (length(red_row) > 0 ) {
  for (row in red_row) {
    points(red_table[row, "bp"], table[i,"segment"], col="red")
  }
}


for (i in c(2:nrow(table))) {
  red_row <- which(red_table$molID == table[i,"molID"])
  if (length(red_row) > 0 ) {
    for (row in red_row) {
      points(red_table[row, "bp"], table[i,"segment"], col="red", pch=20)
    }
  }
  lines(c(table[i, "start"], table[i,"stop"]), rep(table[i,"segment"], 2), type='l', col="blue")
  
}
abline(v=3900000, col="red")
dev.off()


#################################################################################
# Filter out segments with red labels within 10kb of both ends

png("segments_10-50_redlabels_neighbors_10kb.png", height=3000, width=6000, res=300)
first_plot = FALSE
for (i in c(1:nrow(table))) {
  
  red_row <- which(red_table$molID == table[i,"molID"])
  close = FALSE
  
  if (length(red_row) > 0 ) {
    for (row in red_row) {
      red_bp = red_table[row, "bp"]
      seg_start = table[i, "start"]
      seg_end = table[i,"stop"]
      if (! close & ! first_plot) {
        # Check if close
        if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          close = TRUE
          first_plot=TRUE
          plot(c(seg_start, seg_end), rep(table[1,"segment"], 2), type='l', col="blue", xlim=c(0,7500000), ylim=c(0,1944))
          points(red_bp, table[i,"segment"], col="red", pch=20, cex=0.5)
        }
      }
      else if (close) {
          if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          points(red_bp, table[i,"segment"], col="red", pch=20, cex=0.5 )
           }
      }
      else if (! close & first_plot) {
        if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          close = TRUE
          lines(c(seg_start, seg_end), rep(table[i,"segment"], 2), type='l', col="blue")
          points(red_bp, table[i,"segment"], col="red", pch=20, cex=0.5)
      }
      }
    }
    
  }
}
abline(v=3900000, col="red")
dev.off()

#################################################################################
# Filter out segments with red labels within 10kb of both ends
filtered_IDs = c()

png("segments_10-50_redlabels_neighbors_10kb_zoom.png", height=5000, width=4000, res=300)
first_plot = FALSE
for (i in c(1:nrow(table))) {
  
  red_row <- which(red_table$molID == table[i,"molID"])
  close = FALSE
  
  if (length(red_row) > 0 ) {
    for (row in red_row) {
      red_bp = red_table[row, "bp"]
      seg_start = table[i, "start"]
      seg_end = table[i,"stop"]
      if (! close & ! first_plot) {
        # Check if close
        if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          close = TRUE
          first_plot=TRUE
          plot(c(seg_start, seg_end), rep(table[1,"segment"], 2), type='l', col="blue", xlim=c(3600000,4200000), ylim=c(0,1944))
          points(red_bp, table[i,"segment"], col="red", pch=20)
          filtered_IDs <- c(filtered_IDs, table[i, "molID"])
        }
      }
      else if (close) {
        if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          points(red_bp, table[i,"segment"], col="red", pch=20)
        }
      }
      else if (! close & first_plot) {
        if ((abs(seg_start - red_bp) < 10000) | (abs(seg_end - red_bp) < 10000)) {
          close = TRUE
          lines(c(seg_start, seg_end), rep(table[i,"segment"], 2), type='l', col="blue")
          points(red_bp, table[i,"segment"], col="red", pch=20)
          filtered_IDs <- c(filtered_IDs, table[i, "molID"])
        }
      }
    }
    
  }
}
abline(v=3900000, col="red")
dev.off()

print(filtered_IDs)


write(filtered_IDs, "segments_10-50_redlabels_neighbors_10kb_molID_list.txt", sep='\n')






























