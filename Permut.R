#Usage: this program is using permutation test to calculate adjusted p-value.
#Programmer: Binghao Wu
#Time: Oct 10, 2017
#Function #1: calculating t statistics
t_statistics=function (treatment,control) {
  difference_mean=mean(treatment)-mean(control)
  variance_treatment=((var(treatment))^2)/length(treatment)
  variance_control=((var(control))^2)/length(control)
  variance=sqrt(variance_treatment + variance_control)
  t=difference_mean/variance
  t
}
#Function #2: calculating permutated t statistics
Permutation=function (total_value_list,treatment_number) {
  treatment=sample(x=total_value_list,size=treatment_number,replace=FALSE)
  control=total_value_list[-treatment]
  permutation_result=t_statistics(treatment,control)
  permutation_result
}
#Function #3:factorial function
factorial=function (number) { 
  result=1
  for (i in (1:number)) {
    result=result*i
    i=i-1
  }
  result
}
#Function #4: adjusted_pvalue function for single gene
adjusted_p_value_single=function (input,jth_gene,treatment_number,control_number) { 
  #calculate observed t statistics
  matrix_input=as.matrix(input)
  colnames(matrix_input)=NULL
  total_value_list=matrix_input[jth_gene:jth_gene,2:(1+treatment_number+control_number)]
  as.numeric(total_value_list)
  control_group=matrix_input[jth_gene:jth_gene,2:(1+control_number)]  #notice: in the matrix, control must be left,treatment must be right
  treatment_group=matrix_input[jth_gene:jth_gene,(2+control_number):(1+treatment_number+control_number)]
  observed_t_statistics=t_statistics(treatment_group,control_group)
  #calculate permutated t statistics
  n=factorial(treatment_number+control_number)/(factorial(treatment_number))*(factorial(control_number))
  permutated_t_statistics=replicate(n,Permutation(total_value_list,treatment_number))
  #calculate adjusted p value
  B=length(permutated_t_statistics)
  above_observed_number=0
  for(i in 1:B) {
    if (permutated_t_statistics[i]>=observed_t_statistics) {
      above_observed_number=above_observed_number+1
    }
  }
  adjusted_pvalue=above_observed_number/B
  output=matrix(NA,ncol=2)
  colname=c("adjusted_pvalue","difference_mean")
  colnames(output)=colname
  output[,1:1]=adjusted_pvalue
  output[,2:2]=(mean(treatment_group)-mean(control_group))
  output
}

#Example (Notice: in the input data, control must on left, treatment have to on right)
data=read.delim("D:/demo/demo.txt",head=T,sep="\t")
#usage: adjusted_p_value(input,jth_gene,treatment_number,control_number)
adjusted_p_value_single(data,3,3,3)
#Result: 0  0.746895

#calculate genes simultinuously 
#Function #4: adjusted_pvalue function for all genes
adjusted_p_value_total=function (input,treatment_number,control_number) { 
  matrix_input=as.matrix(input)
  colnames(matrix_input)=NULL
  output=matrix(NA,nrow=nrow(matrix_input),ncol=3)
  colname=c("gene_name","adjusted_pvalue","difference_mean")
  colnames(output)=colname
  output[,1:1]=matrix_input[,1:1]
  for(jth_gene in (1:nrow(matrix_input))) {
  #calculate observed t statistics
  total_value_list=matrix_input[jth_gene:jth_gene,2:(1+treatment_number+control_number)]
  as.numeric(total_value_list)
  control_group=matrix_input[jth_gene:jth_gene,2:(1+control_number)]  #notice: in the matrix, control must be left,treatment must be right
  treatment_group=matrix_input[jth_gene:jth_gene,(2+control_number):(1+treatment_number+control_number)]
  observed_t_statistics=t_statistics(treatment_group,control_group)
  #calculate permutated t statistics
  n=factorial(treatment_number+control_number)/(factorial(treatment_number))*(factorial(control_number))
  permutated_t_statistics=replicate(n,Permutation(total_value_list,treatment_number))
  #calculate adjusted p value
  B=length(permutated_t_statistics)
  above_observed_number=0
  for(i in 1:B) {
    if (permutated_t_statistics[i]>=observed_t_statistics) {
      above_observed_number=above_observed_number+1
    }
  }
  adjusted_pvalue=above_observed_number/B
  output[jth_gene:jth_gene,2:2]=adjusted_pvalue
  output[jth_gene:jth_gene,3:3]=(mean(treatment_group)-mean(control_group))
  }
  
  write.csv(output,"D:/demo/output.csv")
  
}

#Example (Notice: in the input data, control must on left, treatment have to on right)
data=read.delim("D:/demo/demo.txt",head=T,sep="\t")
#usage: adjusted_p_value(input,treatment_number,control_number)
adjusted_p_value_total(data,3,3)



