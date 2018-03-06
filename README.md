# Generalized Integrative Principal Component Analysis (GIPCA)

High-dimensional multi-source data are encountered in many fields. Despite recent developments on the integrative dimension reduction of such data, most existing methods cannot easily accommodate data of multiple types (e.g., binary or count-valued). Moreover, multi-source data often have block-wise missing structure, i.e., data in one or more sources may be completely unobserved for a sample. The heterogeneous data types and presence of block-wise missing data pose significant challenges to the integration of multi-source data and further statistical analyses. In this paper, we develop a low-rank method, called Generalized Integrative Principal Component Analysis (GIPCA), for the simultaneous dimension reduction and imputation of multi-source block-wise missing data, where different sources may have different data types. We also devise an adapted BIC criterion for rank estimation.

# Arguments

In `GIPCA` function,<br />
<br />
`data`: A list of two or more linked data matrices . These matrices must have the same row dimension.<br />
<br />
`rankj`: An integer giving the joint rank of the data. <br />
<br />
`ranka`: A vector of integers giving the individual ranks of the data. <br />
<br />
`D`: A vector of integers giving the number of variables in each data set.<br />
<br />
`family`: A vector of characters which sepcify the distributions of each data set in data list.<br />
<br />
`tol`: An integer specifies the tolerance number in the iteration. The default number is 0.1.<br />
<br />
`max.iter`: An integer sepcify the maximum number of iterations. The default number is 500.<br />
<br />
`n.size`: A matrix to specify the number of trials when one of the data set is generated from binomial distribution.<br />
<br />
`lambda`: A penalty term when one of the data set if geenrated from binary distribution. The defual number is 0.<br />
<br />
