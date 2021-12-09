# DLWrap
Deep transfer learning model for disease risk prediction using high-dimensional genomic data.
Building an accurate disease risk prediction model is an essential step towards precision
medicine. While high-dimensional genomic data provides valuable data resources for the
investigations of disease risk, their huge amount of noise and the complex relationships
between predictors and outcomes have brought tremendous analytical challenges. In
this work, we have developed a deep neural network (DNN) based prediction modeling
framework. We first proposed a group-wise feature importance score for feature
selection, where genes harboring genetic variants with both linear and non-linear effects
are efficiently detected. We then developed an explainable transfer-learning based DNN
method, which can directly incorporate information from feature selection and
accurately capture complex predictive effects. The proposed DNN-framework is
computationally efficient and can be applied to genome-wide data. Through extensive
simulation and real data analyses, we have demonstrated that our proposed method can
not only efficiently detect predictive features, but also accurately predict disease risk, as
compared to many existing methods.



## Install
The requirements listed in the requirements.txt are only what's required to install this package so you can use it as a module. They aren't sufficient to actually run all of the scripts and notebooks. In addition, you will need:<br>
```
python=3.6.13
sklearn
numpy==1.19.5
pandas==1.1.5
statistics==1.0.3.5
keras==2.4.3
scipy==1.5.4
tensorflow==2.4.0
```

## Procedure parameters
### Procedure parameters
**TestSelection**:Analysis to be performed. <br>
  * TestSelection=1: Perform screening analysis for one gene; <br>
  * TestSelection=2: Perform screening analysis for multiple genes;<br>
  * TestSelection=3: Perform prediction analysis given screening results; <br>
  * TestSelection=4: Perform both screening and prediction analyses; <br>
  * TestSelection=5: Perform prediction with selected gene.
### Data-related parameters
**train_y_input**: Phenotype of training set: FID, IID, Y and it should have header.**default=None**<br>

**test_y_input**: Phenotype of testing set. It is in the same format as the training set.**default=None**<br>

**binary_outcome**: Whether the outcome is binary: =1 yes binary; =0 no continuous. **Default is continuous (i.e., =0) .**

### Genotypes 
#### 1.only for one gene
**train_x_input**: Genotype of training set,the comma delimited version with the frist 6 columns being FID, IID, PAT, MAT, SEX, PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header. Valid when only one gene is to be screened.**default=None**<br>

**gene_name**: The name of the gene. If not provided, it will be 0. Only valid when screening_index=1.**default=None** <br>

**test_groups**: Groups of tests that need to be performed. If not provided, then all genes will be grouped together. Valid when only one gene is to be screened.**default=None**

#### 2.multiple genes 
**GeneticDataTrain**: The list of training genotype files. Each file contains genotypes of training set, that is the comma delimited version with the frist 6 columns being FID,IID,PAT,MAT,SEX,PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header.**default=None**<br>

**GeneticDataTest**: The list of testing genotype files. Each file contains genotypes of testing set, that is the comma delimited version with the frist 6 columns being FID,IID,PAT,MAT,SEX,PHENOTYPE, and the rest is SNPs. It is in plink .raw format and should have header. The files should be ordered according to the same order in the GeneticDataTrain (i.e., in the same gene order).**default=None**<br>

**geneindexFile**: The list of gene_names, and it should be in the same order as the GeneticDataTrain folder.**default=None**<br>

**test_groupsFile**: The list of groups of test files. It should be in the same order as the GeneticDataTrain folder. Each file contains the grouping for each gene.-----NOT SUPPORTED FOR THE MOMENT!**default=None**

### Output-related parameters
**AssociationDir**: The folder where screening results are saved.**default=None**<br>

**outputPredFile**: The predicted value output.**default=None**

### modules for MLP parameters 
**seed_value**: The random_seed to be set. Default 0. **default=0**<br>

**nunit1**: The number of hidden units in the first hidden layer. **The default is 50. Valid only when screening is performed.**<br>

**nunit2**: The number of hidden units in the second hidden layer. **The default is 10. Valid only when screening is performed.**<br>

**adjusted**: Whether the number of hidden units should be adjusted, **default is True. Valid only when screening is performed.**<br>

**vfold**: The number of cross-validation. **Default is 10. Valid only when screening is performed.**<br>

**nperm**: The number of permutation. **The default is 100. Valid only when screening is performed.**<br>

**reg_tune**: Whether the dropout percentage should be tuned. 1=Yes. 0=No. **The default is 1. Valid only when screening is performed.**<br>

**reg_set**: The dropout percentage. The default is 0.2. If the dropout percentage is to be tuned, then this parameter is ignored. Valid only when screening is performed.<br>
output_level: The level of details of the output. **Valid only when screening is performed.Default is 0.**<br>

 * 0: only useful output is saved. For screening, the p-value for each gene and the model are saved. <br>

 * 1: intermediate level output is saved. For screening, the p-value for each gene, the model, the score are saved. <br>
 
### modules for prediction-related parameters
**alpha**: the significance level for selecting predictive genes. **Valid only prediction is performed.default=0.05.**<br>

**pre_selected_gene_List**: The list of pre-selected genes.**default=None**<br>

**pre_selected_model_list**: The list of models locations for the pre-selected genes.**default=None**

## Usage
* TestSelection=1
  **Perform screening analysis for one gene
```
python ./DLWrap/DLWrapFinal.py  --TestSelection 1 --seed_value 10 --alpha 0.1 --adjusted True --train_y_input YtrainFDG.phen --binary_outcome 0 --AssociationDir Result_List1 --train_x_input GeneTrain_1  --gene_name GeneList_1 --group_index_file group_index

```

* TestSelection=2
  **Perform screening multiple genes 
```
python ./DLWrap/DLWrapFinal.py  --TestSelection 1 --seed_value 10 --alpha 0.1 --adjusted True --train_y_input YtrainFDG.phen --binary_outcome 0 --AssociationDir Result_List1 --GeneticDataTrain GeneTrain_1 --GeneticDataTest GeneTest_1 --geneindexFile GeneList_1

```
* TestSelection=3
  **Prediction only
```
python ./DLWrap/DLWrapFinal.py  --TestSelection 1 --seed_value 10 --alpha 0.1 --adjusted True --train_y_input YtrainFDG.phen --test_y_input YtestFDG.phen --binary_outcome 0 --AssociationDir Result_List1 --GeneticDataTrain GeneTrain_1 --GeneticDataTest GeneTest_1 --geneindexFile GeneList_1 --outputPredFile OutputPred 

```

* TestSelection=4
  **Associaiton + Prediction
```
python ./DLWrap/DLWrapFinal.py  --TestSelection 1 --seed_value 10 --alpha 0.1 --adjusted True --train_y_input YtrainFDG.phen --test_y_input YtestFDG.phen --binary_outcome 0 --AssociationDir Result_List1 --GeneticDataTrain GeneTrain_1 --GeneticDataTest GeneTest_1 --geneindexFile GeneList_1 --outputPredFile OutputPred 

```

## Contributing

PRs accepted.

## License

MIT © Richard McRichface

