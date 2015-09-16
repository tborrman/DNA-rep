# A script to analyze DNA replication origin data by support vector machines

library("e1071");
library("ROCR");

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Read in ori_data_1.6.txt
full_ori_data <- read.table(paste(wkDir, "/ori_data_1.8.txt", sep=""), header=TRUE, sep="\t", comment.char="");

# Classify data as early or late defined by Scott Yang's parameter n from MIM model s.t.
# early = n > median(n)
# late = n <= median(n)

# Get median n value
full_n <- full_ori_data$yang_n
median_n <- median(full_n, na.rm = TRUE);
# Create class vector for labeling early and late origins
class <- sapply(full_n, function(x) {
          if(is.na(x)) {
            return(NA);
          }
          else if(x > median_n) {
            return("early");
          }
          else if(x <= median_n){
            return("late");
          }
});

# Add classification to data
ori_data_class <- cbind(full_ori_data, class);

# Extract all data containing an n parameter
ori_data_clean <- ori_data_class[-which(is.na(full_n)), ];

# Remove remaining origins in rDNA which skew ChIP-seq
ori_data <- ori_data_clean[-which(ori_data_clean$ID == 534),];

# Give 0 to all NAs in rpd3 data
ori_data[which(is.na(ori_data$knott_update_Rpd3_WT_diff)), "knott_update_Rpd3_WT_diff"] <- 0;

# Use 2/3 of data as a training set for svm model
# Remaining 1/3 of data will be the test set


training_total <- round(nrow(ori_data)* (2/3));

training_indices <- sample(1:nrow(ori_data), training_total);
test_indices <- setdiff(1:nrow(ori_data), training_indices);

training_set <- ori_data[training_indices,];
test_set <- ori_data[test_indices,]; 

# Create svm model from training set using just macalpine replicate 2 ChIP seq
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {
  MCM_training <- training_set[c("macalpine_2_MCM_no_mult", "class")];
  MCM_test <- test_set[c("macalpine_2_MCM_no_mult")];
  
  svm_model <- svm(class~., data= MCM_training, kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model);
  
  predict_values <- predict(svm_model, MCM_test, probability = TRUE);

  # Confusion matrix
  confusion_matrix <- table(pred = predict_values, true = test_set$class);

  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late

  TP <- confusion_matrix["early", "early"];
  FP <- confusion_matrix["early", "late"];
  TN <- confusion_matrix["late", "late"];
  FN <- confusion_matrix["late", "early"];

  sensitivity <- TP / (TP + FN);
  specificity <- TN / (TN + FP);

  # Compute ROC curves
  ROC_pred <- prediction(attr(predict_values, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf <- performance(ROC_pred, measure = "tpr", x.measure = "fpr");

  #profilePath <- paste(begPath, "/plot.pdf", sep="");

  #plot(ROC_perf,col="BLUE");

  #dev.off();

  ROC_AUC<-as.numeric(performance(ROC_pred, measure = "auc", x.measure
                                = "cutoff")@ y.values);


  profilePath <- paste(begPath, "/Results/svm/1D/", kern_funct, "_svm_plot_1D.pdf", sep="");
  pdf(profilePath, width=10, height=8);

  #plot(svm_model_2D, MCM_Ku_training);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC, digits = 2), sep = "");
  plot(ROC_perf, col= "BLUE", main = title);

  dev.off();
}

# Create svm model from training set using just ku data
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {
  Ku_training <- training_set[c("donaldson_WT_yku70_diff_plus3", "class")];
  Ku_test <- test_set[c("donaldson_WT_yku70_diff_plus3")];
  
  svm_model_ku_1 <- svm(class~., data= Ku_training, kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_ku_1);
  
  predict_values_ku_1 <- predict(svm_model_ku_1, Ku_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_ku_1 <- table(pred = predict_values_ku_1, true = test_set$class);
  
  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_ku_1["early", "early"];
  FP <- confusion_matrix_ku_1["early", "late"];
  TN <- confusion_matrix_ku_1["late", "late"];
  FN <- confusion_matrix_ku_1["late", "early"];
  
  sensitivity_ku_1 <- TP / (TP + FN);
  specificity_ku_1 <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_ku_1 <- prediction(attr(predict_values_ku_1, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_ku_1 <- performance(ROC_pred_ku_1, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_ku_1 <-as.numeric(performance(ROC_pred_ku_1, measure = "auc", x.measure
                                  = "cutoff")@ y.values);
  
  
  profilePath <- paste(begPath, "/Results/svm/ku_1/", kern_funct, "_svm_plot_1D_ku.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  #plot(svm_model_2D, MCM_Ku_training);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_ku_1, digits = 2), sep = "");
  plot(ROC_perf_ku_1, col= "BLUE", main = title);
  
  dev.off();
}

# Create svm model from training set using just rpd3 data
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {
  Rpd3_training <- training_set[c("knott_update_Rpd3_WT_diff", "class")];
  Rpd3_test <- test_set[c("knott_update_Rpd3_WT_diff")];
  
  svm_model_rpd3_1 <- svm(class~., data= Rpd3_training, kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_rpd3_1);
  
  predict_values_rpd3_1 <- predict(svm_model_rpd3_1, Rpd3_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_rpd3_1 <- table(pred = predict_values_rpd3_1, true = test_set$class);
  
  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_rpd3_1["early", "early"];
  FP <- confusion_matrix_rpd3_1["early", "late"];
  TN <- confusion_matrix_rpd3_1["late", "late"];
  FN <- confusion_matrix_rpd3_1["late", "early"];
  
  sensitivity_rpd3_1 <- TP / (TP + FN);
  specificity_rpd3_1 <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_rpd3_1 <- prediction(attr(predict_values_rpd3_1, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_rpd3_1 <- performance(ROC_pred_rpd3_1, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_rpd3_1 <-as.numeric(performance(ROC_pred_rpd3_1, measure = "auc", x.measure
                                        = "cutoff")@ y.values);
  
  
  profilePath <- paste(begPath, "/Results/svm/rpd3_1/", kern_funct, "_svm_plot_1D_rpd3_1.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  #plot(svm_model_2D, MCM_Ku_training);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_rpd3_1, digits = 2), sep = "");
  plot(ROC_perf_rpd3_1, col= "BLUE", main = title);
  
  dev.off();
}

# Let's add the ku mutant difference in Trep data and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {

  #MCM_col = 82;
  #Ku_col = 63;
  #class_col = 93;

  MCM_Ku_training <- training_set[c("macalpine_2_MCM_no_mult", "donaldson_WT_yku70_diff_plus3", "class")];
  MCM_KU_test <- test_set[c("macalpine_2_MCM_no_mult", "donaldson_WT_yku70_diff_plus3")];

  svm_model_2D <-  svm(class~., data= MCM_Ku_training , kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_2D);

  #plot(svm_model_2D, MCM_Ku_training);

  predict_values_2D <- predict(svm_model_2D, MCM_KU_test, probability = TRUE);

  # Confusion matrix
  confusion_matrix_2D <- table(pred = predict_values_2D, true = test_set$class);

  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late

  TP <- confusion_matrix_2D["early", "early"];
  FP <- confusion_matrix_2D["early", "late"];
  TN <- confusion_matrix_2D["late", "late"];
  FN <- confusion_matrix_2D["late", "early"];

  sensitivity_2D <- TP / (TP + FN);
  specificity_2D <- TN / (TN + FP);

  # Compute ROC curves
  ROC_pred_ku <- prediction(attr(predict_values_2D, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_ku <- performance(ROC_pred_ku, measure = "tpr", x.measure = "fpr");

  #profilePath <- paste(begPath, "/plot.pdf", sep="");

  #plot(ROC_perf,col="BLUE");

  #dev.off();

  ROC_AUC_ku<-as.numeric(performance(ROC_pred_ku, measure = "auc", x.measure
                                      = "cutoff")@ y.values);


  profilePath <- paste(begPath, "/Results/svm/ku/", kern_funct, "_svm_plot_2D_ku.pdf", sep="");
  pdf(profilePath, width=10, height=8);

    plot(svm_model_2D, MCM_Ku_training);
    title = paste(" ROC Plot: AUC = ", round(ROC_AUC_ku, digits = 2), sep = "");
    plot(ROC_perf_ku, col= "BLUE", main = title);

  dev.off();
  
}


# Let's try rpd3 data alone and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {
  
  #MCM_col = 82;
  #Ku_col = 63;
  #class_col = 93;
  
  MCM_rpd3_training <- training_set[c("macalpine_2_MCM_no_mult", "knott_update_Rpd3_WT_diff", "class")];
  MCM_rpd3_test <- test_set[c("macalpine_2_MCM_no_mult", "knott_update_Rpd3_WT_diff")];
  
  svm_model_rpd3 <-  svm(class~., data= MCM_rpd3_training , kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_rpd3);
  
  #plot(svm_model_2D, MCM_Ku_training);
  
  predict_values_rpd3 <- predict(svm_model_rpd3, MCM_rpd3_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_rpd3 <- table(pred = predict_values_rpd3, true = test_set$class);

  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_rpd3["early", "early"];
  FP <- confusion_matrix_rpd3["early", "late"];
  TN <- confusion_matrix_rpd3["late", "late"];
  FN <- confusion_matrix_rpd3["late", "early"];
  
  sensitivity_rpd3 <- TP / (TP + FN);
  specificity_rpd3 <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_rpd3 <- prediction(attr(predict_values_rpd3, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_rpd3 <- performance(ROC_pred_rpd3, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_rpd3 <-as.numeric(performance(ROC_pred_rpd3, measure = "auc", x.measure
                                     = "cutoff")@ y.values);
  
  
  profilePath <- paste(begPath, "/Results/svm/rpd3/", kern_funct, "_svm_plot_2D_rpd3.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  plot(svm_model_rpd3, MCM_rpd3_training);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_rpd3, digits = 2), sep = "");
  plot(ROC_perf_rpd3, col= "BLUE", main = title);
  
  dev.off();
}

# Let's add the Rpd3 dependent data and ku and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");

for (kern_funct in functions) {
  
  #MCM_col = 82;
  #Ku_col = 63;
  #Rpd3_col = 92;
  #class_col = 93;
  
  MCM_Ku_rpd3_training <- training_set[c("macalpine_2_MCM_no_mult", "donaldson_WT_yku70_diff_plus3","knott_update_Rpd3_WT_diff", "class")];
  MCM_KU_rpd3_test <- test_set[c("macalpine_2_MCM_no_mult", "donaldson_WT_yku70_diff_plus3","knott_update_Rpd3_WT_diff")];
  
  svm_model_ku_rpd3 <-  svm(class~., data= MCM_Ku_rpd3_training , kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_ku_rpd3);
  
  #plot(svm_model_2D, MCM_Ku_training);
  
  predict_values_ku_rpd3 <- predict(svm_model_ku_rpd3, MCM_KU_rpd3_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_ku_rpd3 <- table(pred = predict_values_ku_rpd3, true = test_set$class);
  
  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_ku_rpd3["early", "early"];
  FP <- confusion_matrix_ku_rpd3["early", "late"];
  TN <- confusion_matrix_ku_rpd3["late", "late"];
  FN <- confusion_matrix_ku_rpd3["late", "early"];
  
  sensitivity_ku_rpd3 <- TP / (TP + FN);
  specificity_ku_rpd3 <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_ku_rpd3 <- prediction(attr(predict_values_ku_rpd3, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_ku_rpd3 <- performance(ROC_pred_ku_rpd3, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_ku_rpd3 <-as.numeric(performance(ROC_pred_ku_rpd3, measure = "auc", x.measure
                                  = "cutoff")@ y.values);
    
  profilePath <- paste(begPath, "/Results/svm/ku_rpd3/", kern_funct, "_svm_plot_ku_rpd3.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  plot(svm_model_ku_rpd3, MCM_Ku_rpd3_training, macalpine_2_MCM_no_mult ~ knott_update_Rpd3_WT_diff);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_ku_rpd3, digits = 2), sep = "");
  plot(ROC_perf_ku_rpd3, col= "BLUE", main = title);
  
  dev.off();
  
  
}

# Let's add everything we got and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");

for (kern_funct in functions) {
  
  #MCM_col = 82;
  #Ku_col = 63;
  #Rpd3_col = 92;
  #class_col = 93;
  
  MCM_8D_training <- training_set[c("macalpine_2_MCM_no_mult","dsSPD4a_MCM_no_mult", "dsSPD4b_1_MCM_no_mult", "dsSPD8a_MCM_no_mult", "dsSPD8b_MCM_no_mult", "macalpine_1_MCM_no_mult","donaldson_WT_yku70_diff_plus3","knott_update_Rpd3_WT_diff", "class")];
  MCM_8D_test <- test_set[c("macalpine_2_MCM_no_mult", "dsSPD4a_MCM_no_mult", "dsSPD4b_1_MCM_no_mult", "dsSPD8a_MCM_no_mult", "dsSPD8b_MCM_no_mult", "macalpine_1_MCM_no_mult", "donaldson_WT_yku70_diff_plus3","knott_update_Rpd3_WT_diff")];
  
  svm_model_8D <-  svm(class~., data= MCM_8D_training , kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_8D);
  
  #plot(svm_model_2D, MCM_Ku_training);
  
  predict_values_8D <- predict(svm_model_8D, MCM_8D_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_8D <- table(pred = predict_values_8D, true = test_set$class);
  
  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_8D["early", "early"];
  FP <- confusion_matrix_8D["early", "late"];
  TN <- confusion_matrix_8D["late", "late"];
  FN <- confusion_matrix_8D["late", "early"];
  
  sensitivity_8D <- TP / (TP + FN);
  specificity_8D <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_8D <- prediction(attr(predict_values_8D, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_8D <- performance(ROC_pred_8D, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_8D <-as.numeric(performance(ROC_pred_8D, measure = "auc", x.measure
                                           = "cutoff")@ y.values);
  
  profilePath <- paste(begPath, "/Results/svm/8D/", kern_funct, "_svm_plot_8D.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  plot(svm_model_8D, MCM_8D_training, macalpine_2_MCM_no_mult ~ knott_update_Rpd3_WT_diff);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_8D, digits = 2), sep = "");
  plot(ROC_perf_8D, col= "BLUE", main = title);
  
  dev.off();
  
  
}






