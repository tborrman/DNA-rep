# List of functions for data processing
# Notes on merge.with.order:
# Consider merging data X and Y by column C. Denote C column for X and Y as X_c and Y_c respectively.
# Consider merge all.x and keep_order = 1;
# If Y_c is missing a value seen in X_c, the merged data Z will have NAs for missing data from Y
# If X_c is missing a value seen in Y_c, the merged data Z will drop the corresponding data from Y
# If Y_c has two identical values, the merged data Z will have two rows for different Y data corresponding to identical values
# If X_c has two identical values, the merged data z will have corresponding Y data for those identical values
# If 2 identical values for both X_c and Y_c (i.e. 4 rows associated to 1 variable), 
# merged data Z will contain the 4 rows of all merged combinations 
############## function:##################################################################################
# example:     merge.with.order( x.labels, x.vals, by='ref', all.y = T, sort=F ,keep_order = 2) # yay - works as we wanted it to...
# keep_order = 1 keeps x order, keep_order = 2 keeps y order  
merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")
    
    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]		
  }
  
  # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
  # tmp()
  
  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    # if you didn't get "return" by now - issue a warning.
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
  } else {return(merge(x=x,y=y,..., sort = sort))}
}
############################################################################################################################
