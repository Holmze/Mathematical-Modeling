# CH4统计方法建模

## 4.5方差分析

在科学研究中，经常需要分析各种因素对实验指标是否有显著影响，要解决这个问题，一方面是要设计一个实验，使其充分反映各因素的作用，并力求试验次数尽可能少；另一方面要对实验结果的数据进行合理的分析。

### 单因素方差分析

仅考虑一个因素A对于实验指标有无显著影响，可令A取r个值：$A_1,A_2,...,A_r$,称为单因素试验，其结果$x_{ij}$如下：
|序号|1|2|...|$n_i$|
|-|-|-|-|-|
|$A_1$|$x_{11}$|$x_{12}$|...|$x_{1n_i}$|
|$A_2$|$x_{21}$|$x_{22}$|...|$x_{2n_i}$|
|...|...|...|...|...|
|$A_r$|$x_{r1}$|$x_{r2}$|...|$x_{rn_i}$|
假设在$A_i$下的数据$x_{i1},x_{i2},...,x_{in_i}$来自$X_i$~$N(\mu_i,\sigma^2)$，检验如下假设：
$H_0:\mu_1=\mu_2=...=\mu_r,H_1:\mu_1,\mu_2,...,\mu_r不全相等$
检验统计量为：
$F=\frac{\frac{S_A}{r-1}}{\frac{S_e}{n-r}}$~$F(r-1,n-r)$
上式中，
$S_A=\sum_{i=1}^{r}\sum_{j=1}^{n_i}(\overline{x_i}-\overline{x})^2=\sum_{i=1}^{r}(\overline{x_i}-\overline{x})^2$，称为组间差平方和。
$S_e=\sum_{i=1}^{r}\sum_{j=1}^{n_i}(x_{ij}-\overline{x_i})^2$，称为组内差平方和。
其中，$n=\sum_{i=1}^{r}n_i,\\\overline{x_i}=\frac{1}{n}\sum_{j=1}^{n_i}x_{ij},\\\overline{x}=\frac{1}{n}\sum_{i=1}^{r}\sum_{j=1}^{n_i}x_{ij}$
对于给定的显著性水平$\alpha(\alpha=0.0 or 0.05)$，如果$F>F_\alpha(r-1,n-r)$，则拒绝$H_0$，即认为A对实验指标有显著影响。
在实际计算时，可以先对原始数据做如下处理：
$x_{ij}'=\frac{x_{ij}-a}{b}$，不会影响F值的大小。
通常采用如下简便记法：
$T_i=\sum_{j=1}^{n_i}x_{ij} (i=1,2,...,r)$
$T=\sum_{i=1}^{r}\sum_{j=1}^{n_i}x_{ij}$
则有：
$S_A=\sum_{i=1}^{r}\frac{T_i^2}{n_i}-\frac{T^2}{n}$
$S_e=\sum_{i=1}^{r}\sum_{j=1}^{n_i}x_{ij}^2-\sum_{i=1}^{r}\frac{T_i^2}{n_i}$
可得方差分析表：
|方差来源|平方和|自由度|F值|检验p值|
|-|-|-|-|-|
|因素A|$S_A$|r-1|$F=\frac{S_A/(r-1)}{S_e/(n-r)}$||
|误差E|$S_e$|n-r|||
|总和|$S_A+S_e$|n-1|||

### 双因素方差分析

同时考察两个因素A、B对于试验指标有无显著影响，可以让A取r个水平：$A_1,A_2,...,A_r$，让B取r个水平：$B_1,B_2,...,B_r$，在各种水平配合$(A_i,B_j)$下进行实验，称为双因素实验。

#### 1.无交互作用得双因素方差分析

在没中过水平配合$(A_i,B_j)$下做实验，称为无交互作用得双因素实验，检验结果$x_{ij}$如下表：
|因素|$B_1$|$B_2$|...|$B_S$|
|-|-|-|-|-|
|$A_1$|$x_{11}$|$x_{12}$|...|$x_{1s}$|
|$A_2$|$x_{21}$|$x_{22}$|...|$x_{2s}$|
|...|...|...|...|...|
|$A_r$|$x_{r1}$|$x_{r2}$|...|$x_{rs}$|
并假设在水平配合$(A_i,B_j)$下的数据$x_{ij}$来自整体$X_{ij}$~$N(\mu_{ij},\sigma^2),(i=1,2,...,r;j=1,2,...,s)$
检验如下假设：
$H_{0A}:\mu_1.=\mu_2.=...=\mu_r.,H_{1A}:\mu_1.,\mu_2.,...,\mu_r.不全相等$
$H_{0B}:\mu_{.1}=\mu_{.2}=...=\mu_{.r},H_{1B}:\mu_{.1},\mu_{.2},...,\mu_{.r}不全相等$
检验统计量为：
$F_A=\frac{\frac{S_A}{r-1}}{\frac{S_e}{n-r}}$~$F(r-1,n-r)$
$F_B=\frac{\frac{S_B}{r-1}}{\frac{S_e}{n-r}}$~$F(r-1,n-r)$
上式中，
$S_A=\sum_{i=1}^{r}\sum_{j=1}^{s}(\overline{x_i.}-\overline{x})^2=\sum_{i=1}^{r}(\overline{x_i.}-\overline{x})^2$，称为A的组间差平方和。
$S_B=\sum_{i=1}^{r}\sum_{j=1}^{s}(\overline{x_{.j}}-\overline{x})^2=\sum_{i=1}^{s}(\overline{x_{.j}}-\overline{x})^2$，称为B的组间差平方和。
$S_e=\sum_{i=1}^{r}\sum_{j=1}^{n_i}(x_{ij}-\overline{x_i})^2$，称为组内差平方和。
其中，$\overline{x_i.}=\frac{1}{s}\sum_{j=1}^{s}x_{ij},\\\overline{x_{.j}}=\frac{1}{r}\sum_{j=1}^{r}x_{ij},\\\overline{x}=\frac{1}{rs}\sum_{i=1}^{r}\sum_{j=1}^{s}x_{ij}$
对于给定的显著性水平$\alpha(\alpha=0.0 or 0.05)$，
如果$F_A>F_\alpha(r-1,(r-1))(s-1))$，则拒绝$H_{0A}$，即认为A对实验指标有显著影响。
如果$F_B>F_\alpha(r-1,(r-1)(s-1)))$，则拒绝$H_{0B}$，即认为B对时延指标有显著影响。
在实际计算时，可以先对原始数据做如下处理：
$x_{ij}'=\frac{x_{ij}-a}{b}$，不会影响$F_A、F_B$值的大小。
|方差来源|平方和|自由度|F值|检验p值|
|-|-|-|-|-|
|因素A|$S_A$|r-1|$F=\frac{S_A/(r-1)}{S_e/(n-r)}$||
|因素B|$S_B$|r-1|$F=\frac{S_A/(r-1)}{S_e/(n-r)}$||
|误差E|$S_e$|(r-1)(s-1)|||
|总和|$S_A+S_B+S_e$|rs-1|||

#### 2.有交互作用的双因素方差分析

在每一种水平配合$(A_i,B_j)$下重复作$m(m\geq2)$次实验，称之为==有交互的双因素试验==，其结果如下：
![结果](./image/有交互双因素表格.png)
假设在水平配合$(A_i,B_j)$下的数据$x_{ij1},x_{ij2,...,x_{ijm}}$来自总体$X_{ij}$~$N(\mu_{ij},\sigma^2),(i=1,2,...,r;j=1,2,...,s)$。
假设检验如下：
$H_{0A}:\mu_{1.}=\mu_{2.}=...=\mu_{r.},H_{1A}:\mu_{1.},\mu_{2.},...,\mu_{r.}不全相等$
$H_{0B}:\mu_{.1}=\mu_{.2}=...=\mu_{.s},H_{1B}:\mu_{.1},\mu_{.2},...,\mu_{.r}不全相等$
$H_{0AB}:\mu_{ij}全相等，H_{1AB}:\mu_{ij}不全相等$
分别用如下检验统计量：
$F_A=\frac{S_A/(r-1)}{S_e/rs(m-1)}$~$F(r-1,rs(m-1))$
$F_B=\frac{S_B/(r-1)}{S_e/rs(m-1)}$~$F(s-1,rs(m-1))$
$F_{AB}=\frac{S_{AB}/(r-1)(s-1)}{S_e/rs(m-1)}$~$F((r-1)(s-1),rs(m-1))$
其中，
$S_A=\sum_{i=1}^{r}\sum_{j=1}^{s}\sum_{k=1}^{m}(\overline{x_{i.}}-\overline{x})^2=\sum_{i=1}^{r}sm(\overline{x_{i.}}-\overline{x})^2$称为==A的组件差平方和==。
$S_B=\sum_{i=1}^{r}\sum_{j=1}^{s}\sum_{k=1}^{m}(\overline{x_{.j}}-\overline{x})^2=\sum_{i=1}^{r}sm(\overline{x_{.j}}-\overline{x})^2$称为==B的组件差平方和==。
$S_{AB}=\sum_{i=1}^{r}\sum_{j=1}^{s}\sum_{k=1}^{m}(\overline{x_{ij}}-\overline{x_{i.}}-\overline{x_{.j}}+\overline{x})^2=m\sum_{i=1}^{r}\sum_{j=1}^{s}(\overline{x_{ij}}-\overline{x_{i.}}-\overline{x_{.j}}+\overline{x})^2$称为==A*B的组件差平方和==。
$S_e=\sum_{i=1}^{r}\sum_{j=1}^{s}\sum_{k=1}^{m}(x_{ijk}-\overline{x_{ij}})^2$称为==组内差平方和==
这里：
$\overline{x_{i.}}=\frac{1}{sm}\sum_{j=1}^{s}\sum_{k=1}^{m}x_{ijk}$
$\overline{x_{.j}}=\frac{1}{rm}\sum_{i=1}^{r}\sum_{k=1}^{m}x_{ijk}$
$\overline{x_{ij}}=\frac{1}{m}\sum_{k=1}^{m}x_{ijk}$
$\overline{x}=\frac{1}{rsm}\sum_{i=1}^{r}\sum_{j=1}^{s}\sum_{k=1}^{m}x_{ijk}$
对于给定的显著性水平$\alpha(\alpha=0.0 or 0.05)$，
如果$F_A>F_\alpha(r-1,rs(m-1))$，则拒绝$H_{0A}$，即认为A对实验指标有显著影响。
如果$F_B>F_\alpha(r-1,rs(m-1))$，则拒绝$H_{0B}$，即认为B对时延指标有显著影响。
如果$F_{AB}>F_\alpha((r-1)(s-1),rs(m-1))$，则拒绝$H_{0AB}$，即认为A与B的交互效应对实验指标有显著影响。

在实际计算时，可以先对原始数据做如下处理：
$x_{ij}'=\frac{x_{ij}-a}{b}$，不会影响$F_A、F_B、F_{AB}$值的大小。

**方差分析表**：
|方差来源|平方和|自由度|均方|检验F值|p|
|-|-|-|-|-|-|
|因素A|$S_A$|r-1|$\frac{S_A}{r-1}$|$F=\frac{S_A/(r-1)}{S_e/rs(m-1)}$||
|因素B|$S_B$|r-1|$\frac{S_B}{s-1}$|$F=\frac{S_B/(s-1)}{S_e/rs(m-1)}$||
|交互AB|$S_{AB}$|(r-1)(s-1)|$\frac{S_{AB}}{(r-1)(s-1)}$|$F=\frac{S_A/(r-1)(s-1)}{S_e/rs(m-1)}$||
|误差E|$S_e$|rs(m-1)|$\frac{S_e}{rs(m-1)}$||
|总和|$S_A+S_B+S_{AB}+S_e$|rsm-1|||

## 4.6回归分析

### 一、多元线性回归

如果随机变量Y与固定变量$x_1,x_2,...,x_m$之间有明显的线性关系，即：$Y=b_0+b_1x_1+b_2x_2+...+b_mx_m+\epsilon,\epsilon$~$N(0,\sigma^2)$称为m元线性回归模型。

#### 1.模型中的参数估计——回归方程的建立

设通过试验/历史资料得到的观测数据$(x_{i1},x_{i2},...,x_{im},y_i),(i=1,2,...,n)$
令
$$
Y=
\begin{pmatrix}
    y_1\\
    y_2\\
    .\\
    .\\
    .\\
    y_n
\end{pmatrix}
$$
$$
X=
\begin{pmatrix}
    1&x_{11}&x_{12}&...&x_{1m}\\
    1&x_{21}&x_{22}&...&x_{2m}\\
    .&.&.&...&.\\
    .&.&.&...&.\\
    .&.&.&...&.\\
    1&x_{n1}&x_{n2}&...&x_{nm}\\
\end{pmatrix}
$$
$$B=
\begin{pmatrix}
    b_0\\
    b_1\\
    .\\
    .\\
    .\\
    b_m
\end{pmatrix}
$$
由最下二乘估计，可得：$\hat{B}=(X^TX)^{-1}X^TY$
称$\hat{Y}=\hat{b_0}+\hat{b_1}x_{i1}+\hat{b_2}x_{i2}+...+\hat{b_m}x_{im}(i=1,2,...,n)$

#### 2.显著性检验

##### (1)回归模型的显著性检验

即检验假设$H_0:b_1=b_2=...=b_m=0,H_1:b_i不全为0$，令$S_R=\sum_{i=1}^{n}(\hat{y_i}-\overline{y})^2,S_e=\sum_{i=1}^{n}(y_i-\hat{y_i})^2$
检验统计量：
$F=\frac{S_R/m}{S_e/(n-m-1)}$~$F(m,n-m-1)$对于一个小概率$\alpha$，若$F>F_\alpha(m,n-m-1)$,则拒绝$H_0$，认为所见的回归方程正确。