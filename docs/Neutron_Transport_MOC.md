# 中子输运方程特征线法简述

稳态的中子输运方程可写为：
\[
    \Omega \cdot \nabla\phi + \Sigma_t \phi = Q
\]
其中$Q$为总源项，$Q=Q_s+Q_f+S$ 。
$Q_s$为散射源，$Q_f$为裂变源，$S$为外源。
考虑方程在沿中子运动方向的特征线上的表达形式，中子输运方程可写为：
\[
    \frac {d\phi}{ds} + \Sigma_t \phi = Q
\]
其中$s$为沿特征线方向的距离。
此时方程化为以$s$为自变量的一阶线性非齐次常微分方程。
设$Q$已知，则可求得其通解为：
\[
    \phi = e^{-\Sigma_t s}(\int_0^s Q(s') e^{-\Sigma_t s'} ds' + C)
\]
对能量、位置、角度进行离散，设源项各向同性，设一网格内的各源项为常数。
则第$g$能群，第$j$方向的第$k$条特征线段上的中子输运方程可写为：
\[
    \frac {d\phi_{g,j,k}}{ds} + \Sigma_{t,g} \phi_{g,j,k} = Q_{g,k}
\]
设入射处$\phi^{in}$已知，则可求得其解为：
\[
    \phi_{g,j,k} (s)=\phi_{g,j,k}^{in} e^{-s\Sigma_{t,g}} + \frac {Q_{g,k}}{\Sigma_{t,g}} (1-e^{-s\Sigma_{t,g} })
\]
设该特征线段的长度为$l$，则可得出射处：
\[
    \phi_{g,j,k}^{out} = \phi_{g,j,k}^{in} e^{-l\Sigma_{t,g}} + \frac {Q_{g,k}}{\Sigma_{t,g}} (1-e^{-l\Sigma_{t,g}})
\]
该特征线段的平均中子角通量密度为:
\[
    \overline\phi_{g,j,k} = \frac{Q_{g,k}}{\Sigma_{t,g}} + \frac{(\phi_{g,j,k}^{in} - \frac{Q_{g,k}}{\Sigma_{t,g}})(1-e^{-l\Sigma_{t,g}})}{\Sigma_{t,g}l}
\]
出射处的值可作为下一特征线段的入射值。
对各能群、方向、特征线段上的中子输运方程进行计算，则可得中子角通量密度$\phi(r,\Omega,E)$，进而可算得中子通量密度$\Phi(r,E)$和有效增殖因子$k_{eff}$。