# 動学BLP
以下のモデルは[Dubé,Fox,Su(2009)](https://www.nber.org/papers/w14991)で導出されている動学BLPモデルです。 動学BLPではタイプ $r$ の消費者が、 $t$ 期に財 $j \in J$ を購入する意思決定をモデル化したものです。消費者は先読み行動(Forward Looking Behavior)をとるため、財がカメラやゲーム機の場合、消費者は価格が下がるのを待つために購入タイミングをずらします。動学BLPではそのような消費者の意思決定を扱います。さらに消費者をセグメント $r \in R$ で分けることで消費者がセグメントごとで異なる価格感応度や選好をもつことを表現できます。消費者の買い控えする時の効用はベルマン方程式を用いて以下の様に表されます。

## モデル

$$
V_0^r(p; \theta^r)= 
$$

$$
\delta \int_{ \epsilon,p } \\{ \max \\{ V_0^r(p_{j,t+1}; \theta^r ) + \epsilon_{0,t}, \max_{j} \\{ \beta_j^r - \alpha^r p_{j,t+1 } + \xi_{j,t} + \epsilon_{j,t} \\}  \\}  \\} d_{\epsilon}F(\epsilon)d_{p,\xi}F(p,\xi)
$$

$$
=\delta \int_{p,\xi} log\left(exp(V_0^r(p_{j,t+1}; \theta^r))+
{\large\Sigma_{j}} exp(\beta_{j}^r-\alpha^r p_{j,t+1} + \xi_{j,t}\right)d_{p,\xi}F(p,\xi)
$$

$p$ は価格、 $\theta^r$ はタイプ $r$ 消費者のパラメータ、 $\beta,\alpha$ はそれぞれ財への選好と価格感応度です。 $\epsilon$ は攪乱項でタイプI極値分布に従います。価格は以下の様に遷移します。

$$
p_{j,t+1} = ( \rho_{j,0} \space \rho_{j,1} \space \rho_{j,2} \space \cdots) 
\begin{pmatrix}
1 \\
p_{j,t} \\
p_{j,t}^2 \\
\vdots \\
\end{pmatrix} + \psi_{i,t} = p_{t}^{'} \rho_{j}+\psi_{j,t}\hspace{5pt}
for\hspace{5pt}j=1,2,...J
$$

すると上のベルマン方程式は以下の様になります。

$$
V_0^r(p; \theta^r)= 
$$

$$
=\delta \int_{\psi,\xi} log\left(exp(V_0^r( p_{t}^{'} \rho_{j}+\psi_{j,t}; \theta^r))+
{\large\Sigma_{j}} exp(\beta_{j}^r
-\alpha^r ( p_{t}^{'} \rho_{j}+\psi_{j,t}) + \xi_{j,t}
\right)d_{\psi,\xi}F(\psi,\xi)
$$

財 $j \in J$ のシェア計算は以下の通りです。 $\lambda_{t,r}$ は消費者のタイプの割合です。タイプによってパラメータ $\beta,\alpha$ が異なります。

$$
s_j(p_t;\theta) = 
\sum_{r=1}^{R}\lambda_{t,r} \frac{\exp(\beta_{j}^r-\alpha^r p_{j,t}+\xi_{j,t})}
{ \exp( v_0^r(p_{j,t}; \theta^r ) ) + \sum_{k=1 \space to \space J}\exp( \beta_k^r - \alpha^r p_{k,t} + \xi_{k,t} ) }, \space j = 1,2,...,J
$$

マーケットサイズ推移とタイプ割合推移

$$
M_t^r =
\begin{cases} 
M^r \lambda_{r} & ,t = 1 \\ 
M_{t-1}^r s^0(p_{t-1},\xi_{t-1};\theta^r) &, t > 1
\end{cases}
$$

$$
\lambda_{t,r} = 
\begin{cases}
\lambda_r &, t=1 \\
\frac{\lambda_r}{M_t^r} &, t>1
\end{cases}
$$

## 推定

$Y_t = (p_t,S_t)$ とすると、

$$
f_Y (Y_t; \theta, \rho, \Omega) = \frac{1}{(2\pi)^J |\Omega|^{1/2}} \exp\left(-\frac{1}{2} u_t'\Omega^{-1}u_t\right)|J_{t,u→Y}|
$$

where

$$
u_t =
\begin{bmatrix}
\psi_{t} \\
\xi_{t}
\end{bmatrix}=
\begin{bmatrix}
p_t-p_{t-1}^{'}\rho \\
s^{-1}(p_t,S_t;\theta)
\end{bmatrix}
$$

$$
\max_{\{\theta,\rho,\Omega\}} \prod_{t=1}^T f_Y (Y_t; \theta, \rho, \Omega)
$$

ヤコビ行列は以下の様になります

$$
J_{t,u \rightarrow Y} = \begin{bmatrix} \frac{\partial \psi_t}{\partial p_t} & \frac{\partial \psi_t}{\partial S_t} \\ \frac{\partial \xi_t}{\partial p_t} & \frac{\partial \xi_t}{\partial S_t} \end{bmatrix}
$$

$$
G(S_t,\xi_t:\theta) = s(p,\xi_t; \theta) - S_t=0
$$

$$
J_{t,\xi \rightarrow S} = -\begin{bmatrix}\frac{\partial G}{\partial \xi_{t}}
\end{bmatrix}^{-1}
\begin{bmatrix}
\frac{\delta G}{\delta S_{t}}\\
\end{bmatrix}=
\begin{bmatrix}
\frac{\delta s}{\delta \xi_{t}}\\
\end{bmatrix}^{-1}
$$

$$
\begin{equation}
\frac{\delta s_{j,t}} {\delta \xi_{k,t}} = 
\begin{cases} 
{ \Sigma }_{r} \lambda _{r,t} s_j ( p_t, \xi ; \theta^r ) ( 1-s_j ( p _t, \xi_t; \theta^r )) \hspace{5pt} ,\text{if} \space j=k\\
{ \Sigma } _{r} \lambda _{r,t} s_j ( p_t, \xi_t; \theta^r ) s_k ( p, \xi_t; \theta^r ) \hspace{5pt} ,\text{otherwise} \\
\end{cases}
\end{equation}
$$

### MPEC推定

$$
\begin{equation}
\begin{aligned}
\max_{\theta,\rho,\Omega,\xi,v} \prod_{t=1}^{T} \frac{1}{(2\pi)^{\frac{J}{2}}|\Omega|^\frac{1}{2}} \exp \left(-\frac{1}{2}u_t'\Omega^{-1}u_t \right) |J_{t,u \rightarrow Y} | 
\end{aligned}
\end{equation}
$$

$\text{subject to}$

$$
\begin{equation}
\begin{aligned}
s(\xi_t;\theta) = S_t, \space  \space \forall t = 1,...,T \\
\end{aligned}
\end{equation}
$$

$\text{and}$

$$
v_0^r(p) \\
= \delta log \left( \exp(v_0^{r}(p^{'}\rho_j+\psi)) + ... 
\sum_j^{}\exp(\beta_j^r -\alpha^r (p^{'}\rho_j +\psi) + \xi_{j,t})\right)dF_{\psi,\xi}(\psi, \xi), \\
\forall d ∈ D, r = 1,...,R.
$$


### NFXP推定
内ループ(価値関数反復)

$$ 
T(V_0^r) - V_0^r = 0 \space \text{for} \space r = 1,2,..,R
$$

$T$ はベルマンオペレータ。

中ループ(縮小写像)

$$
\xi_{\text{new}} = \xi_{\text{old}} + log \space S - log \space s(\xi_{\text{old}};\theta) 
$$

外ループ（最適化）

$$
\begin{equation}
\begin{aligned}
\max_{\theta,\rho,\Omega,\xi,v} \prod_{t=1}^{T} \frac{1}{(2\pi)^{\frac{J}{2}}|\Omega|^\frac{1}{2}} \exp \left(-\frac{1}{2}u_t'\Omega^{-1}u_t \right) |J_{t,u \rightarrow Y} | 
\end{aligned}
\end{equation}
$$
