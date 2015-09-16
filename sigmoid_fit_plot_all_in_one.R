# A file to fit a sigmoid to Shankar's Nanostring data and compute Trep and goodness of fit
# Note all plots are put into one pdf differing from original sigmoid_fit.R which plots separate pdfs

# Set paths
begPath <- "/Users/User/Research/DNArep";
# Open data table
data <-read.table(paste(begPath, "/Data/shankar/Dec5&18-1.txt", sep=""), header=TRUE, sep="\t", comment.char="");

# Add time header to dataset initialize variables
colnames(data)[1] <- "time";
origin_name <- c();
Trep <- c();
residuals_sum <- c();
colors <- c();

# Extract time column
minutes <- data$time;

# Sigmoid function with adjustable parameters
sigmoid <- function(params, minutes) {
  params[1] + (params[4]/(1 + exp(-params[2] * (minutes-params[3]))));
  
}
# Plot result
profileFile <-"sigmoidfit_Dec5&18-1.pdf";
profilePath <- paste(begPath,"/Data/shankar/vis/", profileFile, sep="");
pdf(profilePath, width=6, height=6);
  title <- "Sigmoidal fits for nanostring data";
  # Fit a sigmoid to each origin
  origin <- 2;
  ori_name <- colnames(data)[origin];
  copy_number <- data[,origin];
  # Get non-linear least square estimates for parameters of sigmoid function
  fitModel = nls(copy_number ~  a+(d/(1 + exp(-b*(minutes-c)))), start=list(a=1,b=.5, c=40, d=1), trace=TRUE, control = nls.control(warnOnly = TRUE));
  params = coef(fitModel);
  # Make list of minute values for function plotting
  new_minutes <- seq(0,80);
  # Get new model fit values
  fit_values = sigmoid(params, new_minutes);
  plot(new_minutes, fit_values, type= 'l', ylim=c(0.8, 2.3), pch = 20, main = title, ylab = "Copy number", xlab = "minutes", col = origin);
  # Save color numbers for legend (colors are numbered in R)
  colors <- c(colors, origin);
  # Save origin name
  origin_name <- c(origin_name, ori_name);
  Trep <- c(Trep, params[3]);
  # Get sum of residuals
  residuals_sum <- c(residuals_sum, sum(abs(residuals(fitModel)))); 

  for (origin in c(3:ncol(data))) {
    ori_name <- colnames(data)[origin];
    copy_number <- data[,origin];
    # Get non-linear least square estimates for parameters of sigmoid function
    fitModel = nls(copy_number ~  a+(d/(1 + exp(-b*(minutes-c)))), start=list(a=1,b=.5, c=40, d=1), trace=TRUE, control = nls.control(warnOnly = TRUE));
    params = coef(fitModel);
    # Make list of minute values for function plotting
    new_minutes <- seq(0,80);
    # Get new model fit values
    fit_values = sigmoid(params, new_minutes);
    lines(new_minutes, fit_values, pch=20, col=origin);
    # Save color numbers for legend (colors are numbered in R)
    colors <- c(colors, origin);
    # Save origin name
    origin_name <- c(origin_name, ori_name);
    # Get Trep
    Trep <- c(Trep, params[3]);
    # Get sum of residuals
    residuals_sum <- c(residuals_sum, sum(abs(residuals(fitModel))));  
  }
  legend("topleft", origin_name, col=colors, lty=1, bty = "n");
  dev.off();

# Fit a sigmoid to each origin
#for (origin in c(2:ncol(data))) {
  #ori_name <- colnames(data)[origin];
  #copy_number <- data[,origin];
  # Get non-linear least square estimates for parameters of sigmoid function
  #fitModel = nls(copy_number ~  a+(d/(1 + exp(-b*(minutes-c)))), start=list(a=1,b=.5, c=40, d=1), trace=TRUE, control = nls.control(warnOnly = TRUE));
  #params = coef(fitModel);
  # Make list of minute values for function plotting
  #new_minutes <- seq(0,80);
  # Get new model fit values
  #fit_values = sigmoid(params, new_minutes);
  # Plot result
  #profileFile <- paste("sigmoidfit_", ori_name, ".pdf", sep="");
  #profilePath <- paste(begPath,"/Data/shankar/vis/", profileFile, sep="");
  #pdf(profilePath, width=6, height=6);
   # title <- paste("Copy number (blue) and sigmoidal fit (red) for ", ori_name, sep="");
    #plot(minutes, copy_number, main = title,  ylab= "Copy number", pch = 20, col = "blue");
    #lines(new_minutes, fit_values, pch=20, col="red");
  #dev.off();
  # Save origin name
  #origin_name <- c(origin_name, ori_name);
  # Get Trep
  #Trep <- c(Trep, params[3]);
  # Get sum of residuals
  #residuals_sum <- c(residuals_sum, sum(abs(residuals(fitModel))));  
#}

# Output Trep and residual table for all origins
output_table <-  cbind(origin_name, Trep, residuals_sum);
write.table(output_table, paste(begPath,"/Data/shankar/Trep_Residual_matrix.txt", sep=""), row.names=F, sep='\t');
