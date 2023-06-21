# Radiomic_Macrophage

### 中文说明

数据处理主要使用jupyter notebook(python版本3.6~3.9均可)和R语言；绘图只使用R语言。

提交的代码涉及三个主要部分，每个部分均有对应的jupyter数据处理和R语言绘图代码：

- 影像组学预处理【0~2】
- 机器学习建模、验证 【4~7】
- GSVA、cibersort、estimate【8~9】

另外，数据概览和三线图绘图的R代码也提供了【3】。每一部分代码均给出了案例数据供代码运行测试，其中JSPH数据集均匿名化处理。

注：template-GliomaROI.mrb是MRI预处理阶段的3Dslicer模板文件，里面有用于重采样的标准模板数据，供参考和下载使用。

### English description

Data processing mainly uses Jupyter Notebook (Python versions 3.6-3.9 can be used) and [R version 4.2.3](https://www.r-project.org/) .The drawing only uses the R.

The submitted code involves three main parts, each with corresponding jupyter notebook and R codes for data processing and drawing:

- MRI preprocessing, feature extraction and feature cleaning【0~3】
- machine learning model construction and validation【4~7】
- GSVA、cibersort、estimate【8~9】

In addition, the R code for data overview and three line diagram drawing also provides【3】.Each part of the code provides demo data for the code running test, in which the JSPH dataset is anonymized.

