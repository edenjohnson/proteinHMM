import os

# create the HMM with 80% of the basic data set sequences
basic_command = "java -jar phmm2.jar basic_training_set.fa"
os.system(basic_command)

# create the HMM with 80% of the related data set sequences
related_command = "java -jar phmm2.jar related_training_set.fa"
os.system(related_command)