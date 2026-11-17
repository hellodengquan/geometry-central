# 向量热方法中的切向量传输：局部基、方向与强度的语义分析

> 分析对象：geometry-central（分支 `task/50`）中 `VectorHeatMethodSolver` 及其依赖的联络算子。
> 主要源码：
> - `include/geometrycentral/surface/vector_heat_method.h`
> - `src/surface/vector_heat_method.cpp`（下称 **VHM.cpp**）
> - `src/surface/intrinsic_geometry_interface.cpp`（下称 **IGI.cpp**）
> - `include/geometrycentral/utilities/vector2.ipp`（下称 **V2.ipp**）
> - 对照实现：`src/pointcloud/point_cloud_heat_solver.cpp`、`src/surface/polygon_mesh_heat_solver.cpp`

---

## 0. 结论速览

1. **复数坐标是"带基坐标"，不是几何向量本身。** 每个顶点有自己独立的二维切平面基（`halfedgeVectorsInVertex`，IGI.cpp:367–391），跨顶点的数值差异同时包含"真实方向变化"和"基的旋转"。同一几何向量换一组合法局部基，复坐标只相差一个单位复数旋转（规范变换）；联络 Laplacian 的构造保证整个求解是**规范协变**的，解码回几何向量后不变。
2. **方向与强度在算法中被人为拆成两个通道。** 方向通道先把每个源向量归一化为单位复数，经复联络热方程扩散，再逐点归一化；强度通道走独立的标量热方程插值，最后两者相乘（VHM.cpp:274、308、314–335）。
3. **归一化只保留相位（方向），热扩散本身不保长。** 联络的跨面耦合项是纯单位旋转（精确平行传输保长），但 `(M + tL_conn)^{-1}` 是衰减平滑算子；输出长度是后处理"重新施加"的，不是传输保持的。
4. **退化情形没有防护。** 曲面版对零向量和相消点直接调用无保护的 `Vector2::normalize()`（V2.ipp:72–75），零除产生 `NaN`/`Inf` 并污染整个顶点场；点云版改用了带阈值的 `normalizeCutoff`（point_cloud_heat_solver.cpp:144），说明这是已知缺口而非理论必然。
5. **复数坐标数值相等 ≠ 方向相同。** 跨顶点比较方向必须先用 `transportVectorsAlongHalfedge` 做基变换；多源场在相消区的归一化方向由数值残差主导，没有几何意义。

---

## 1. 跨面局部基与联络算子的构造

### 1.1 顶点切平面的复坐标约定

`computeHalfedgeVectorsInVertex()`（IGI.cpp:367–391）为每个顶点建立一套**内在的**二维坐标系：

```cpp
// IGI.cpp:378-389
double coordSum = 0.0;
Halfedge firstHe = v.halfedge();
Halfedge currHe = firstHe;
do {
  halfedgeVectorsInVertex[currHe] = Vector2::fromAngle(coordSum) * edgeLengths[currHe.edge()];
  coordSum += cornerScaledAngles[currHe.corner()];
  ...
} while (currHe != firstHe);
```

要点：

- 角度原点锚定在该顶点的"第一条出半边" `v.halfedge()`，逆时针排列；
- 相邻扇区的夹角用 `cornerScaledAngles`（IGI.cpp:148–163），即把真实角和归一化到内点 `2π`、边界点 `π`。这是把顶点邻域**展开（unfolding）到平面**的角坐标，与 3D 嵌入无关——`IntrinsicGeometryInterface` 下不存在三维基；
- 因而每个顶点的基是**各自独立选取**的，两个顶点的复坐标 `z_i, z_j` 不能直接比较。
- 该量要求流形网格（implicit twin），否则直接抛异常（IGI.cpp:369–371）。

面内另有一套基 `halfedgeVectorsInFace`（IGI.cpp:291–336）：每个三角形单独放到平面上，首边沿实轴。

### 1.2 沿半边的平行传输：单位复数旋转

`computeTransportVectorsAlongHalfedge()`（IGI.cpp:397–414）：

```cpp
// IGI.cpp:407-412
Vector2 vecA = halfedgeVectorsInVertex[heA];   // 边在尾顶点基中的方向（背离尾点）
Vector2 vecB = halfedgeVectorsInVertex[heB];   // twin：在头顶点基中（背离头点，几何反向）
Vector2 rot = unit(-vecB / vecA);              // 复数除法 + 强制单位模
transportVectorsAlongHalfedge[heA] = rot;      // 尾 → 头
transportVectorsAlongHalfedge[heB]  = rot.inv();
```

`Vector2` 的乘除就是复数乘除（V2.ipp:8–13）。几何含义：沿边方向的向量 `vecA` 平行传输到另一端后应对准 `-vecB`，故传输系数 `rot = (-vecB)/vecA`，再经 `unit()` 保证 `|rot|=1`。**联络本身是纯旋转，精确算术下严格保长。** 面间版本 `transportVectorsAcrossHalfedge`（IGI.cpp:341–359）构造方式相同，只是改用面基。

### 1.3 顶点联络 Laplacian

`computeVertexConnectionLaplacian()`（IGI.cpp:525–545）：

```cpp
// IGI.cpp:532-541
Vector2 rot = transportVectorsAlongHalfedge[he.twin()]; // 头 j → 尾 i 的旋转
triplets.emplace_back(iTail, iTail,  weight);
triplets.emplace_back(iTail, iTip,  -weight * rot);
```

即

$$
(L_{\mathrm{conn}} z)_i = \sum_{j\sim i} w_{ij}\bigl(z_i - \mathrm{rot}_{j\to i}\, z_j\bigr),$$

先把邻点的值旋转进本点基再作差。非对角元是复数、对角元是实数，矩阵 Hermitian（Delaunay 条件下正定，VHM.cpp:57–80 据此选择 Cholesky 或通用方阵求解器）。

### 1.4 规范变换（换基）分析——输出应如何解释

设顶点 $v$ 换一组合法基，新旧基相差旋转角 $\alpha_v$，同一几何向量的坐标变换为

$$z'_v = g_v z_v,\qquad |g_v|=1\quad(g_v=e^{-i\alpha_v}).$$

由 1.2 的几何定义，传输系数同步变换为

$$\mathrm{rot}'_{j\to i} = g_i\,\mathrm{rot}_{j\to i}\,g_j^{-1}.$$

代入联络 Laplacian：

$$(L'_{\mathrm{conn}} z')_i
 = \sum_j w_{ij}\bigl(g_i z_i - g_i\,\mathrm{rot}\,g_j^{-1}g_j z_j\bigr)
 = g_i (L_{\mathrm{conn}} z)_i.$$

即 $L_{\mathrm{conn}}$ **规范协变**。质量矩阵 $M$ 是实对角阵，与 $g_v$ 可交换，因此热算子 $M+tL_{\mathrm{conn}}$ 保持同一协变性：若源项按 $b'_v=g_v b_v$ 换基，则解满足 $z'_v=g_v z_v$。

**判断（问题一）：同一几何向量改用另一组合法局部坐标描述时——**

- 输入坐标必须按新基旋转后传入；求解输出的每个复坐标也整体乘上该顶点的规范因子；
- 用各自顶点基解码出的几何切向量完全相同，换基不改变任何几何结果；
- 因此输出的正确解释方式是"**顶点 $v$ 切平面内、以该顶点 `halfedgeVectorsInVertex`  convention 为基的方向分量**"，而不是一组全局坐标系下的二维数；
- 任何两个顶点上的复数值**不能直接判等**。要判断方向是否一致，必须先做传输对齐，例如比较 $z_i$ 与 $\mathrm{rot}_{j\to i}z_j$（即把 $j$ 的值转进 $i$ 的基），再比幅角。`computeLogMap_VectorHeat` 中 `logDir = radialSol/horizontalSol`（VHM.cpp:545）正是用两个单位复数量的商消去规范因子、提取相对角度；`computeLogMap_AffineLocal` 中 `logCoord / parallel_section[v]`（VHM.cpp:378）是同一思想的显式基变换。

---

## 2. 向量热传输的完整数据流

入口 `transportTangentVectors(sources)`（VHM.cpp:243–340）。

### 2.1 右端项：方向归一化 + 重心散布

```cpp
// VHM.cpp:271-303
std::complex<double> unitVec = Vector2::fromComplex(vec).normalize(); // ① 逐源归一化
magnitudeSources.emplace_back(point, vec.norm());                     // ② 强度另行登记
SurfacePoint facePoint = point.inSomeFace();                          // ③ 转到某个相邻面
...
dirRHS[vInd] += w * unitVec;                                          // ④ 重心权重 w 散布
```

- 每个源在进入方向通道前被压成单位复数，原始模长存入 `magnitudeSources`；
- `SurfacePoint`（顶点/边/面内点）经 `inSomeFace()`（surface_point.ipp:39–86）转为面内重心坐标，三个权重即 `faceCoords.x/y/z`，顶点源退化为 `(1,0,0)`；
- 同一顶点收到多个源贡献时用 `+=` 累加（VHM.cpp:285 等），**相消可以在 RHS 阶段就发生**。

### 2.2 复热方程求解

```cpp
// VHM.cpp:53, 308
SparseMatrix<std::complex<double>> vectorOp = massMat.cast<...>() + shortTime * Lconn;
Vector<std::complex<double>> vecSolution = vectorHeatSolver->solve(dirRHS);
```

`shortTime = tCoef * meanEdgeLength²`（VHM.cpp:16–21）。扩散把源方向沿曲面平滑铺开，跨面旋转由 $L_{\mathrm{conn}}$ 的单位复系数承担。注意扩散解的**模长携带距离衰减与多源混合信息，不保留源长度**。

### 2.3 单源：归一化后乘回原长度

```cpp
// VHM.cpp:314-323
double targetNorm = std::get<1>(sources[0]).norm();
vecSolution = (vecSolution.array() / vecSolution.array().abs()) * targetNorm;
```

逐点除以自身模长得到单位方向，再统一乘源向量模长。长度恒定是**构造性强制**的，而非热方程保持的。

### 2.4 多源：方向归一化 × 标量强度插值

```cpp
// VHM.cpp:328-334
VertexData<double> interpMags = extendScalar(magnitudeSources);
for (Vertex v : mesh.vertices()) {
  Vector2 dir = Vector2::fromComplex(vecSolution[...]).normalize();
  result[v] = dir * interpMags[v];
}
```

方向通道与强度通道彻底分离：方向来自复场的相位，强度来自独立的标量最近源扩展。

### 2.5 标量强度插值的内部结构

`extendScalar`（VHM.cpp:171–227）对同一组源解**两次**标量热方程：

```cpp
// VHM.cpp:216-221
Vector<double> dataSol      = scalarHeatSolver->solve(dataRHS);      // 权重 × 强度值
Vector<double> indicatorSol = scalarHeatSolver->solve(indicatorRHS); // 仅权重
Vector<double> interpResult = dataSol.array() / indicatorSol.array();
```

比值消去热核的局部密度差异，得到以"最近源 Voronoi 扩展 + 短时光滑"形式的插值。由于 `indicatorRHS` 与强度值无关（只累加重心权重，VHM.cpp:194 等），在含源的连通分量上热核严格为正，**这个除法本身对零强度是稳健的**——零强度源使 `dataSol` 趋于 0，但分母不受影响。

---

## 3. 专题判断

### 3.1 多源方向相消或强度为零时，归一化是否仍有意义？退化如何传播？

**归一化在相消区没有几何意义，且现有曲面实现无任何防护。** 传播链条如下：

**(a) 零向量源在入口即产生 NaN。** `Vector2::normalize()` 无零保护：

```cpp
// V2.ipp:72-75
inline Vector2 Vector2::normalize() const {
  double r = 1. / std::sqrt(x*x + y*y);   // 零向量 → r = Inf
  return *this * r;                        // 0 * Inf = NaN
}
```

VHM.cpp:274 对每个源调用它。源向量为零 → `unitVec` 为 NaN → `dirRHS` 带 NaN → 线性求解结果整列 NaN → 输出全场 NaN。注意此时强度通道登记的 `vec.norm()=0`（VHM.cpp:277）本身合法，问题完全出在方向通道。

**(b) 扩散场相消点触发零除。** 两种情形：

- 多个源在某顶点 RHS 处复数累加恰好（或近似）相消（VHM.cpp:285 的 `+=`）；
- 网格含**多个连通分量**：无源分量上 `dirRHS=0`，齐次解 `vecSolution=0`，单源路径 VHM.cpp:318 的 `array()/array().abs()` 与多源路径 VHM.cpp:332 的 `normalize()` 都执行 `0/0 → NaN`。

**(c) 近相消比精确相消更隐蔽。** `|vecSolution|` 很小时，其辐角由各源贡献的数值残差决定，对网格精度、`tCoef`、求解器误差极度敏感；归一化把这点微小残差放大成一个"堂而皇之"的单位方向。它不平行于任何一个源，也不具有拓扑定义，纯属数值假象。

**(d) 零强度不能消除方向垃圾。** 多源最终结果是 `dir * interpMags[v]`：若插值强度严格为 0，`0 × NaN/Inf = NaN`，垃圾照样传播；若强度只是很小，得到一个模长很小但方向错误的向量，调用方若只检查非零就会误用其方向。

**(e) 对照证据：点云版有防护、曲面版没有。** 点云求解器使用带阈值的版本：

```cpp
// point_cloud_heat_solver.cpp:142-145
Vector2 val{...};
dirInterp(i) = normalizeCutoff(val);   // 模长 ≤ 阈值时按长度 1 处理，不产生 NaN
```

`normalizeCutoff`（V2.ipp:80–85）在 `len <= mag` 时退回单位长 1。多边形网格版（polygon_mesh_heat_solver.cpp:165、170）与曲面版一样是无保护的裸 `normalize()`。这说明退化处理是**实现选择上的缺口**，不是算法要求如此。

**结论：** 强度为零是可正常插值的标量事实；但"零强度处方向无定义"这一语义在现有实现中没有表达载体——既不返回零向量也不标记无效，而是依赖浮点异常值传播。调用方必须自行把 `Vector2::isFinite()`/`isDefined()`（V2.ipp:157–160）与模长阈值检查作为验收条件。

### 3.2 平行传输方向经归一化后保留了什么信息？

只保留**相位（辐角）信息**，即在每个顶点自身切基中的方向角：

- 被丢弃的是复扩散解的模长——它混合了距离衰减、热核密度、多源权重，与"被传输向量的长度"无关；
- 保留下来的辐角在短时间极限 $t\to 0$ 下趋近真正的平行传输方向（VHM 论文的核心结论），有限 $t$ 下是平滑近似；
- 规范协变性（1.4 节）保证这个辐角随局部基同步旋转，即它编码的是**几何方向**而非坐标假象；
- 对数图进一步展示了归一化场的正确用法：径向场与水平场各自归一化后相除（VHM.cpp:545），用商的辐角取相对方向，模长信息完全改由 Poisson 解出的距离 `distance[v]` 提供（VHM.cpp:534–546）。

### 3.3 标量强度插值如何参与最终结果？能否共同视为保长的几何传输？

**单源情形：** 输出在每个有限顶点上满足 `|result[v]| = targetNorm`（VHM.cpp:316–318），长度处处等于源长度。但这是**重新投影**出来的常数长度，热扩散过程本身不保长；方向是平行传输的短时间近似。因此严格表述应是"**近似平行传输的方向场 + 精确重施加的常数长度**"，在结果有限、网格连通的前提下可实用地称为保长传输，但它不是由传输算子保持的保长，且在零源/不连通网格上不成立。

**多源情形：** 不存在一个被整体传输的几何向量。最终场

$$\text{result}[v] = \underbrace{\frac{z(v)}{|z(v)|}}_{\text{复热场相位：混合方向}} \times \underbrace{\frac{(M+tL)^{-1}b_{\text{mag}}}{(M+tL)^{-1}b_{\text{ind}}}}_{\text{标量最近源插值强度}},$$

是两个独立扩展算子的乘积：在每个源的 Voronoi 影响区内趋近该源的传输方向与强度，在影响区边界两侧通道同时平滑过渡。它是一个**连续的向量扩展/插值场**，不是任何单个向量的平行传输，"保长"无从谈起——长度本来就应随位置切换到不同源的强度。文档措辞 "extend transport the vectors ... nearest source point"（docs/.../vector_heat_method.md:66）也印证这是扩展语义。

**共通限制：** 两个通道都带热方法固有的短时光滑畸变（文档同节自陈 "slightly smoothed-out field"），归一化不消除方向畸变，只消除幅度畸变。

---

## 4. 调用方输入约定

| 约定 | 证据与说明 |
|---|---|
| 流形三角网格、implicit twin | `computeHalfedgeVectorsInVertex` 对一般 `SurfaceMesh` 抛异常（IGI.cpp:369–371）；头文件亦有 TODO 注明应限定 manifold（vector_heat_method.h:64–65） |
| `Vector2` 是源点切平面内的**复坐标**，不是 3D 向量 | 内在几何无三维基；需 3D 表达时要像 polygon 版那样另配 `vertexTangentBasis` 解码（polygon_mesh_heat_solver.cpp:196） |
| 顶点源的基方向 | 实轴沿 `v.halfedge()`，角度按 `cornerScaledAngles` 逆时针计（IGI.cpp:382–389）；调用方若持有另一套合法基，须先按规范因子旋转输入 |
| 边/面内源的坐标语义 | 经 `inSomeFace()` 以重心权重散布到面的三个顶点（VHM.cpp:279–302），依赖各顶点角展开基在同一平面三角形内的一致性；输入向量按所在面展开平面的切坐标解释 |
| 空源列表 | 标量返回全 `NaN`，向量返回 `Vector2::undefined()`（VHM.cpp:172–174、245–247），调用方需显式拦截 |
| 零向量源 | 无防护，产生全场 NaN（见 3.1）；调用方应在传入前过滤或改传显式零并绕过方向求解 |
| 时间参数 | `shortTime = tCoef * meanEdge²`，构造后不可改（VHM.cpp:16–21，头文件:51–52） |
| 非 Delaunay 网格 | cotan 权出现 `< -1e-6` 时退化为通用方阵求解器（VHM.cpp:57–80），符号差网格上方向场质量无保证 |
| 多连通分量 | 无源分量输出为 NaN（零除），需按分量分别求解或事后掩码 |
| 输出容器 | `VertexData<Vector2>`，值在**各顶点自己的基**中；跨顶点比较前必须经 `transportVectorsAlongHalfedge` 对齐 |

---

## 5. 结果判据：现有保证 vs 需补充的处理

### 5.1 源码可支撑的现有保证

| 性质 | 保证来源 |
|---|---|
| 跨面传输是纯旋转、精确保长的联络 | `unit(-vecB/vecA)` 强制单位模（IGI.cpp:409），逆传输取共轭倒数（IGI.cpp:412，V2.ipp:124） |
| 求解与输出规范协变（换基不变的几何语义） | 联络 Laplacian 的 `-w·rot` 构造（IGI.cpp:538–540）与 1.4 节推导；log map 的商/除基操作（VHM.cpp:378、545）是工程印证 |
| 单源输出长度恒定（在有限值处） | 显式重投影 `z/|z| * targetNorm`（VHM.cpp:318） |
| 标量强度插值对零强度稳健 | 指标函数 RHS 与值无关、热核正性保证分母非零（VHM.cpp:181、194、221） |
| 输出场连续、短时光滑 | 椭圆算子 `M + tL` 的解（VHM.cpp:39、53） |
| 空输入有确定（虽无效）输出 | NaN/undefined 哨兵（VHM.cpp:173、246） |

### 5.2 源码不保证、需调用方或补丁补充的处理

| 缺口 | 现场行为 | 建议判据/处理 |
|---|---|---|
| 零向量源 | 入口 `normalize()` 零除 → 全场 NaN（VHM.cpp:274 + V2.ipp:72–75） | 传入前判 `vec.norm()==0`：强度为零的源应只进标量通道并令方向输出置零，或改用 `normalizeCutoff` 风格防护 |
| 扩散相消 / 近相消 | `z/|z|` 零除或放大残差方向（VHM.cpp:318、332） | 归一化前设模长阈值（相对局部质量/RHS 尺度），低于阈值标记方向未定义；可直接复用 `normalizeCutoff`（V2.ipp:80–85） |
| 不连通分量 | 无源分量零解除得 NaN | 按连通分量掩码，或验收时逐点 `isFinite()`/`isDefined()`（V2.ipp:157–160） |
| 零强度区的方向垃圾 | `0 × NaN = NaN`，小强度携带错误方向（VHM.cpp:333） | 强度低于阈值时结果显式置零，不读取方向 |
| 跨顶点方向比较 | 复数值直接相等只代表"同基坐标相等"，基不同则无意义 | 先 `transportVectorsAlongHalfedge[he.twin()] * z_j` 对齐到同一基，再比 `arg()`（V2.ipp:126）；**严禁以复数坐标数值相等判定方向相同** |
| 多源场的语义 | 实为最近源扩展，非全局平行传输 | 验收按影响分区进行：源邻域内比对单源传输精度，边界区只要求连续性；不要对全场施加"长度守恒"判据 |
| 有限时间近似误差 | 方向带热光滑偏差，`tCoef` 敏感 | 用已知方向（如 log map 径向场）做角度残差检查；需要高精度时切换 `LogMapStrategy::AffineAdaptive` 或内蕴 Delaunay 三角化（文档 vector_heat_method.md:105–129） |

### 5.3 最小验收流程建议

1. **输入侧**：过滤空源与零向量源；确认网格流形且已压缩；记录每个源所用切基。
2. **求解侧**：逐点检查 `result[v].isFinite()`，对 NaN/Inf 顶点按连通分量或相消区掩码隔离。
3. **方向校验**：任选边两端顶点，把一端结果用 `transportVectorsAlongHalfedge` 转入另一端基后比较幅角差；单源情形该差值应仅反映热方法光滑误差。
4. **长度校验**：单源时 `|result[v]|` 应在数值误差内恒等于源模长（不满足即触发了退化路径）；多源时改为校验强度通道与 `extendScalar` 输出一致。
5. **基无关性校验**：人为旋转任一顶点基的输入约定并相应旋转输出，解码后的几何向量应不变——这是规范协变性的直接测试，也是区分"坐标变化"与"真实方向变化"的根本判据。

---

## 附：对数图中的同类模式（佐证设计意图）

- `computeLogMap_VectorHeat`：径向、水平两个复场分别归一化（VHM.cpp:495、509），源点手动置零（VHM.cpp:496），最终以单位复数之商取方向、Poisson 距离取长度（VHM.cpp:545–546）——与传输接口"方向/强度双通道"同构；
- `computeLogMap_AffineLocal`：`logCoord / parallel_section[v]`（VHM.cpp:378）显式把各顶点坐标换回源点基，说明库内消费方向场时从不裸比复数，而是先消规范；
- `computeLogMap` 对边/面源的混合（VHM.cpp:651–683）同样先算 `rot` 基变换再线性组合，可作为调用方处理跨基数据的范式。
