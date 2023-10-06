# UCP
All optimization problems are implemented in Python using scikit-learn API and in Matlab using the YALMIP toolbox with Gurobi solver 8.1.1. We provide two case studies including IEEE 6-bus and 118-bus test systems, each residing in its designated folder. The description of files and data are given below.

### Overview of files
#### IEEE 6-bus test system:

- "Data_6bus.mat"
  - This includes the technical data for conventional units, transmission lines, wind and load profiles. You need this data file to run the unit commitment problem.

- "AC_SOCP_6bus"
  - This file runs the unit commitment problem for each sample which is denoted by xx_new and the optimal commitment of generating units are represented by y. The unit commitment problem is modeled as mixed-integer second-order cone program.   

- "classification_SVM_crossvalidation_i"
  - After running the unit commitment problem for each sample and getting labels y, this file runs linear support vector machine to find classifiers. The regularization parameters are tuned by cross validation with 4 folds.

- "classification_SVM_crossvalidation_i_k"
  - This file extends the previous file to kernelized support vector machine.

- "main.py"
  - This file is in Python which corresponds to the files named "classification_SVM_crossvalidation_i" and "classification_SVM_crossvalidation_i_k". We use this files to visualize the performance of classifiers.
  
- "UC.py"
  - We use this file to score the performance of classifiers based on the hing loss function. You can do the same with the files named "classification_SVM_crossvalidation_i" and "classification_SVM_crossvalidation_i_k". Here two different test datasets are provided, each with 1, 000 samples, one is more similar to the training dataset named "baised_test_dataset", while the other one named "unbaised_test_dataset" is less similar. 


#### IEEE 118-bus test system: 
We have files similar to the IEEE 6-bus test system. In this folder, load and wind profiles and the features used for training and testing are provided in seperate files.
