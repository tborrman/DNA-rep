# Trying to learn support vector machine algorithms from package 'e1071'
# Simple examples to apply to replication origins


library("e1071");

# EXAMPLE 1 ---------------------------------------------------------------

# CLASSIFY 1 DIMENSIONAL DATA
# MCM ChIP SEQ and Scott Yang's "n" parameter for 20 origins
# n parameter classified origins into early and late such that early = n > median(n) and late = n <= median(n)
# Following standard svm protocol early = -1 and late = 1

# CLASSIFY NEW UNKNOWN DATA 
# Create linearly separable data
early_MCM_ChIP <- sample(101:200, 10);
late_MCM_ChIP <- sample(1:100, 10);

MCM_ChIP <- c(early_MCM_ChIP, late_MCM_ChIP);
n <- c(rep(-1, 10), rep(1,10));
EX_1_data <- data.frame(MCM_ChIP, n);

# Create svm object
model <- svm(EX_1_data[1],EX_1_data[2], cost = 1, gamma = 1, type = "C-classification");
summary(model);

# Create 10 unknown origins and classify as early or late
new_unknown <- sample(1:200, 10);
pred_new_unknown <- predict(model, new_unknown);

# CLASSIFY SUBSET OF DATA
subset <- sample(MCM_ChIP, 5);
pred_subset <- predict(model, subset);


# EXAMPLE 2 ---------------------------------------------------------------

# CLASSIFY 1 DIMENSIONAL DATA
# CLASSIFY NEW UNKNOWN DATA
# This time use gene expression as example to classify cancer patients
# -1 type is gene expression near 0 and normal patient
# 1 type is low or high gene expression and indicative of cancer patient

# Create linearly inseparable data
EX_2_high_exp <- sample(5:30, 5);
EX_2_low_exp <- sample(-30:-5, 5);
EX_2_norm <- sample(-4:4, 10, replace = TRUE);


EX_2_exp <- c(EX_2_low_exp, EX_2_norm, EX_2_high_exp);
EX_2_type<- c(rep(1, 5), rep(-1,10), rep(1,5));
EX_2_data <- data.frame(EX_2_exp, EX_2_type);

# Create svm object
EX_2_model <- svm(EX_2_data[1],EX_2_data[2], cost = 1, gamma = 1, type = "C-classification");
summary(EX_2_model);

# Create 10 unknown origins and classify as early or late
EX_2_unknown <- sample(-30:30, 10, replace=TRUE);
EX_2_pred_unknown <- predict(EX_2_model, EX_2_unknown);



# EXAMPLE 3 ---------------------------------------------------------------

# CLASSIFY 1 DIMENSIONAL DATA
# cLASSIFY NEW UNKNOWN DATA
# MCM and n inseparable case
# Construct same setup as in example 1, but this time noise will be introduced such that
# some of the origins will have a randomly assigned label 

# How much noise? Adjust this parameter to introduce random labels for origins with ChIP seq data
EX_3_noise <- 40; # Started to notice mistakes when 30 randomly assigned origins were added to original 20 origin training set

# CLASSIFY NEW UNKNOWN DATA
EX_3_early_MCM_ChIP <- sample(101:200, 10);
EX_3_late_MCM_ChIP <- sample(1:100, 10);
EX_3_random_MCM_ChIP <- sample(1:200, EX_3_noise);

EX_3_MCM_ChIP <- c(EX_3_early_MCM_ChIP, EX_3_late_MCM_ChIP, EX_3_random_MCM_ChIP);
EX_3_n <- c(rep(-1, 10), rep(1,10), sample(c(-1,1), EX_3_noise, replace = TRUE));
EX_3_data <- data.frame(EX_3_MCM_ChIP, EX_3_n);

# Create svm object
# Note high cost = strict margin
# low cost = soft margin (default = 1)

EX_3_model <- svm(EX_3_data[1],EX_3_data[2], cost = 10, gamma = 1, type = "C-classification");
summary(EX_3_model);

# Create 10 unknown origins and classify as early or late
EX_3_unknown <- sample(1:200, 10);
EX_3_pred_unknown <- predict(EX_3_model, EX_3_unknown);







