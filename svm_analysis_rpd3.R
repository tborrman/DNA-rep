# A script to analyze DNA replication origin data by support vector machines

library("e1071");
library("ROCR");

# Set paths
begPath <- "/Users/User/Research/DNArep";
wkDir <- paste(begPath, "/Data", sep="");

# Read in ori_data_1.6.txt
full_ori_data <- read.table(paste(wkDir, "/ori_data_1.6.txt", sep=""), header=TRUE, sep="\t", comment.char="");

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

# Use 2/3 of data as a training set for svm model
# Remaining 1/3 of data will be the test set

training_total <- round(nrow(ori_data)* (2/3));

training_indices <- sample(1:nrow(ori_data), training_total);
test_indices <- setdiff(1:nrow(ori_data), training_indices);

training_set <- ori_data[training_indices,];
test_set <- ori_data[test_indices,]; 

# Create svm model from training set

svm_model <- svm(training_set$dsSPD4_MCM_count, training_set$class , cost = 1, gamma = 1, type = "C-classification");
summary(svm_model);


predict_values<- predict(svm_model, test_set$dsSPD4_MCM_count);

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


# Let's add the ku mutant difference in Trep data and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");
for (kern_funct in functions) {

  MCM_col = 82;
  Ku_col = 63;
  class_col = 93;

  MCM_Ku_training <- training_set[c(MCM_col, Ku_col, class_col)];
  MCM_KU_test <- test_set[c(MCM_col, Ku_col)];

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
  ROC_pred <- prediction(attr(predict_values_2D, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf <- performance(ROC_pred, measure = "tpr", x.measure = "fpr");

  #profilePath <- paste(begPath, "/plot.pdf", sep="");

  #plot(ROC_perf,col="BLUE");

  #dev.off();

  ROC_AUC<-as.numeric(performance(ROC_pred, measure = "auc", x.measure
                                      = "cutoff")@ y.values);


  profilePath <- paste(begPath, "/Results/svm/", kern_funct, "_svm_plot_2D.pdf", sep="");
  pdf(profilePath, width=10, height=8);

    plot(svm_model_2D, MCM_Ku_training);
    title = paste(" ROC Plot: AUC = ", round(ROC_AUC, digits = 2), sep = "");
    plot(ROC_perf, col= "BLUE", main = title);

  dev.off();
  
  
}

# Let's add the Rpd3 dependent data and train again!
functions <- c("linear", "polynomial", "radial", "sigmoid");

for (kern_funct in functions) {
  
  MCM_col = 82;
  Ku_col = 63;
  Rpd3_col = 92;
  class_col = 93;
  
  MCM_Ku_training <- training_set[c(MCM_col, Ku_col, Rpd3_col, class_col)];
  MCM_KU_test <- test_set[c(MCM_col, Ku_col, Rpd3_col)];
  
  svm_model_3D <-  svm(class~., data= MCM_Ku_training , kernel = kern_funct, cost = 1, type = "C-classification", probability = TRUE);
  summary(svm_model_3D);
  
  #plot(svm_model_2D, MCM_Ku_training);
  
  predict_values_3D <- predict(svm_model_3D, MCM_KU_test, probability = TRUE);
  
  # Confusion matrix
  confusion_matrix_3D <- table(pred = predict_values_3D, true = test_set$class);
  
  # Compute sensititivity = TP/ (TP + FN) and specificity = TN / (TN + FP)
  # such that Positive = early and Negative = late
  
  TP <- confusion_matrix_3D["early", "early"];
  FP <- confusion_matrix_3D["early", "late"];
  TN <- confusion_matrix_3D["late", "late"];
  FN <- confusion_matrix_3D["late", "early"];
  
  sensitivity_3D <- TP / (TP + FN);
  specificity_3D <- TN / (TN + FP);
  
  # Compute ROC curves
  ROC_pred_3D <- prediction(attr(predict_values_3D, "probabilities") [,"early"], test_set$class == "early" );
  ROC_perf_3D <- performance(ROC_pred_3D, measure = "tpr", x.measure = "fpr");
  
  #profilePath <- paste(begPath, "/plot.pdf", sep="");
  
  #plot(ROC_perf,col="BLUE");
  
  #dev.off();
  
  ROC_AUC_3D<-as.numeric(performance(ROC_pred_3D, measure = "auc", x.measure
                                  = "cutoff")@ y.values);
    
  profilePath <- paste(begPath, "/Results/svm/", kern_funct, "_svm_plot_3D.pdf", sep="");
  pdf(profilePath, width=10, height=8);
  
  plot(svm_model_3D, MCM_Ku_training, dsSPD4_MCM_count ~ knott_sig_diff);
  title = paste(" ROC Plot: AUC = ", round(ROC_AUC_3D, digits = 2), sep = "");
  plot(ROC_perf, col= "BLUE", main = title);
  
  dev.off();
  
  
}






