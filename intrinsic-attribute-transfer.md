# 内蕴翻边改变三角划分后，分片线性属性的迁移方案评估

分析对象：geometry-central，分支 `task/35`，commit `019669d`（"Improve marching triangles"）。
本文件中的每条结论都标注来源：**[代码]** 由源码结构或代数恒等式直接导出；**[实测]** 在本机编译该 commit（Clang, Release, 无 SuiteSparse → Eigen `SparseLU`）后用真实 API 复现，数值已记录在第 6 节。

记号：

- `A` = `inputMesh`（原三角划分），`B` = `intrinsicMesh`（内蕴三角划分），`CS` = 公共细分。
- `P_A: R^{V_A} → R^{V_CS}`、`P_B: R^{V_B} → R^{V_CS}`：`CommonSubdivision::interpolationMatrixA/B`。
- `M`：CS 上的 `vertexGalerkinMassMatrix`（`transfer_functions.h:29-31`）。
- "迁移"指把 A 上的 PL 标量属性 `f` 变成 B 上的节点值 `g`（或反向）。

---

## 0. 结论摘要

1. **迁移的几何对应由数据记录给出，而不是由测地追踪给出[代码]。** CS 的每个顶点是一条 `CommonSubdivisionPoint{intersectionType, posA, posB, orientation}`（`common_subdivision.h:31-37`），同时携带该点在两侧曲面上的位置；`P_A`/`P_B`、`interpolateAcrossA/B`、`copyFromA/B` 与 `AttributeTransfer` 的全部实现都只读这条记录（`common_subdivision.cpp:100-182, 316-344`；`common_subdivision.ipp:8-46`；`transfer_functions.cpp:57-105`）。`equivalentPointOnIntrinsic/OnInput` 是另一条独立路径，迁移不使用它。
2. **但"记录质量"依赖后端，且 signpost 在开边界网格上会引入退化单元[实测]。** 对"在开边界网格的内边中点插点"这一操作，signpost 后端把新顶点在 A 上的位置记成 **面内点**（几何上恰在边上，偏差 3.3e-16），并额外产生一个几乎重合的 `EDGE_TRANSVERSE` 点，于是 CS 变成 6 顶点/6 面、其中 **2 个面面积为 0**；整数坐标后端记成 **边上点**（`tEdge = 0.5`），CS 为 5 顶点/4 面、无退化面。两者迁移出的数值仍一致（差异 ≤ 3.3e-16），但下游若按 CS 面遍历/可视化，signpost 的结果含零面积单元。闭合网格上不发生（E12）。
3. **Pointwise 与 L2 的保证互补且互斥[代码+实测]。** A→B→A 方向：Pointwise 逐位恢复 A 的全部顶点值（实测精确相等）；L2 则在 B 上给出不同的函数，往返不恒等（实测 `[0,1,0,1] → [0.25,0.5,0.25,0.5]`）。两者只当原函数落在两侧 PL 空间的交集时才同时满足。
4. **L2 的面积守恒是代数恒等式，与对应精度无关[代码+实测]。** 因为 `P_A 1 = P_B 1 = 1`，`1ᵀMP_A = 1ᵀMP_B·T`（实测逐列成立）。它守恒的是"以所用 M 为测度的积分"。
5. **M 本身可能不是真实面积[实测]。** 整数坐标后端默认做 mollification（`mollifyEPS = 1e-5`）。在 h=1e-3 的针状三角形上，`AttributeTransfer(intTri)` 用的 M 给出面积 `9.13e-4`，而真实面积是 `5.0e-4`（**1.83 倍**）。同时 `intTri.faceAreas` 在构造后仍报告 **mollification 之前**的 `5.0000000001e-4`——缓存未失效，内部自相矛盾。
6. **`AttributeTransfer` 的生命周期没有保护[代码+实测]。** 它持有 CS 的引用并把 `P_A`/`P_B`/`M` 拷贝出来；一旦迁移变异三角划分，CS 被 `reset()` 释放。实测：先做一次 L2 迁移、再翻一条边、再迁移 ⇒ **段错误**；Pointwise 路径则静默返回按旧 CS 尺寸计算的结果。
7. **退化/歧义情形的边界明确[代码]，且多数无测试[实测/检索]。** `orientation` 在 signpost 后端恒为 `false`；共享边配对靠启发式搜索；`constructMeshData` 遇退化面只往 stderr 打一行并跳过；`recoverTraceTValues` 把累计长度归一化到 `[0,1]`，在被 `maxIters` 截断的追踪上静默给出错误位置。
8. **向量/方向场不在这两个 API 的保证范围内[代码]。** `SurfacePoint::interpolate` 逐分量线性插值，不做平行移动；A、B 顶点切标架在翻边后不保证一致；`AttributeTransfer` 只支持 `double`。

---

## 1. 现有方案与数据流

| 能力 | 入口 | 位置 | 输入 → 输出 | 依赖 CS |
| --- | --- | --- | --- | --- |
| 逐点采样 A→B | `IntrinsicTriangulation::sampleFromInput` | `intrinsic_triangulation.ipp:22-29` | `VertexData<T>` on A → on B | 否（只读 `vertexLocations`） |
| 逐点限制 B→A | `IntrinsicTriangulation::restrictToInput` | `intrinsic_triangulation.ipp:31-40` | on B → on A | 否 |
| CS 上点采样 | `CommonSubdivision::interpolateAcrossA/B` | `common_subdivision.ipp:8-26` | on A/B → on CS | 是 |
| CS 面常量拷贝 | `CommonSubdivision::copyFromA/B` | `common_subdivision.ipp:28-46` | `FaceData<T>` on A/B → on CS | 是 |
| 插值矩阵 | `interpolationMatrixA/B` | `common_subdivision.cpp:100-182` | `P_A`, `P_B` | 是（需 `constructMesh`） |
| 逐点迁移 | `AttributeTransfer::transferAtoB_Pointwise` 等 | `transfer_functions.cpp:57-70, 82-95` | `VertexData<double>` | 是 |
| L2 最优迁移 | `AttributeTransfer::transferAtoB_L2` 等 | `transfer_functions.cpp:72-80, 97-105` | `VertexData<double>` | 是 |

两条路径的数值关系：

- **[代码]** `transferAtoB_Pointwise` 取 `result[p.posB.vertex] = p.posA.interpolate(valuesOnA)`（`transfer_functions.cpp:63-69`）；`sampleFromInput` 取 `vertexLocations[v].interpolate(...)`（`intrinsic_triangulation.ipp:25-27`）。CS 中 B 顶点的 `posA` 就是 `vertexLocations[vB]`（`signpost_intrinsic_triangulation.cpp:777-806`；`integer_coordinates_intrinsic_triangulation.cpp:199-224`），故两者在 B 的所有顶点上给出同一数值，差别只是前者要构建 CS。
- **[实测]** 在翻边后的正方形上，`sampleFromInput` 与 `transferAtoB(..., Pointwise)` 都输出 `[0,1,0,1]`，`restrictToInput` 输出 `[0,1,0,1]`，与该结构判断一致。

---

## 2. 问题一：顶点、边如何对应到同一曲面位置；保留了哪些连续性

### 2.1 顶点

CS 顶点分五类（`common_subdivision.h:23-29`）[代码]：

| `CSIntersectionType` | 含义 | `posA` | `posB` | 迁移中的作用 |
| --- | --- | --- | --- | --- |
| `VERTEX_VERTEX` | A、B 公共顶点 | `SurfacePoint(vA)` | `SurfacePoint(vB)` | Pointwise 唯一"写值"的点 |
| `EDGE_VERTEX` | B 的新顶点落在 A 的边上 | `SurfacePoint(eA, t)` | `SurfacePoint(vB)` | signpost 会额外登记到 `pointsAlongA`（`signpost_intrinsic_triangulation.cpp:791-798`）；整数坐标不登记（`integer_coordinates_intrinsic_triangulation.cpp:213-216`） |
| `FACE_VERTEX` | B 的新顶点落在 A 的面内 | `SurfacePoint(fA, bary)` | `SurfacePoint(vB)` | 按重心坐标插值 |
| `EDGE_TRANSVERSE` | B 的边横穿 A 的边 | `SurfacePoint(eA, t)` | `SurfacePoint(eB, t)` | 两侧 `pointsAlong` 都记录 |
| `EDGE_PARALLEL` | 两边的"平行"伪交点（HACK） | `SurfacePoint(eA, 0.5)` | `SurfacePoint(eB, 0.5)` | 记录共享边；**不是**真实交点 |

共享顶点按**索引**对应：signpost 用 `getIndex()`（`signpost_intrinsic_triangulation.cpp:107, 120, 133, 152-157`），整数坐标显式假设"A 的顶点在 B 中同索引"（`integer_coordinates_intrinsic_triangulation.cpp:134-135, 199-223`）。这就是文档前置条件 "Every vertex of `meshA` must also be a vertex of `meshB`"（`docs/.../common_subdivision.md:11`）[代码]。

### 2.2 "同一曲面位置"的实际质量：后端相关[实测]

对"在 A、B 共享的内边中点插点"这一完全相同的操作，两个后端给出的记录不同：

| 后端 | 新顶点的 `posA` | 几何位置误差 | CS |
| --- | --- | --- | --- |
| `SignpostIntrinsicTriangulation` | **`Face`** 型（面点） | 3.3e-16（到对角线） | 6 顶点 / 6 面，**2 个面面积 = 0** |
| `IntegerCoordinatesIntrinsicTriangulation` | `Edge` 型，`tEdge = 0.5` | 0 | 5 顶点 / 4 面，无退化面 |

机制（E4 的完整 CS dump）：signpost 的 `resolveNewVertex` 沿 A 面内的方向**追踪**来定位（`signpost_intrinsic_triangulation.cpp:726-737`），追踪的终点是"面内点"（`trace_geodesic.h:30-35` 说明默认不在边上结束），于是即便几何上恰好落在共享边上，记录仍是 `FACE_VERTEX`，其 `posA.faceCoords = (0.5, 0.5, 3.3e-16)`。随后把 B 的内边沿 A 追溯时，这条内边因该点不在边上而有微小偏移，额外登记了一个 `EDGE_TRANSVERSE`，其 `posB.tEdge ≈ 3.5e-16`（几乎正好落在同一个几何点）。两个几何重合但索引不同的 CS 顶点 ⇒ **2 个零面积面**（E4 中面积 `0, 0`，其余 4 个各 `0.25`）。整数坐标后端对共享边（`normalCoordinates < 0`）走 `splitInteriorEdge` 专门分支，直接返回 `SurfacePoint(inputEdge, tInsertion)`（`integer_coordinates_intrinsic_triangulation.cpp:1450-1471`），没有这个重复点。

后果：**[实测 E4/探针18]** 两者插值出的数值一致 —— 插点情形下 PW 与 L2 在两个后端上差异 ≤ 3.3e-16（`SurfacePoint::interpolate` 只看几何位置，类型只决定分支）；CS 总面积仍为 1，且**没有** `"degree-2 face"` 警告。但 CS 的**单元结构不同**：signpost 版本含 2 个零面积面，任何"按 CS 面遍历/累加/可视化"的下游都会看到退化单元，而库不作任何提示。

**稳定性与适用范围[实测 E12]**：该现象与 `tSplit` 取值无关（`t = 0.25, 0.3, 0.5, 0.7` 均产生 6V6F 与 2 个零面积面），也不是"共享边"独有——只要拆分点位于开边界网格的一条内边上就会触发。在**闭合**网格（四面体）上先 `flipToDelaunay()` 再 `splitEdge`，CS 为 6V8F、**无**零面积面。仓库测试全部使用闭合网格（`fox.ply`、`cat_head.obj`）且不检查零面积面，因此这条缺陷不会被现有测试发现。

### 2.3 边

- **[代码]** `pointsAlongA[eA]` / `pointsAlongB[eB]` 是沿边有序的 CS 点表（含端点），排序键为 `posA.tEdge`（`signpost_intrinsic_triangulation.cpp:855-870`）。
- **[代码]** `intersectionsA/B` 返回去掉端点后的交点数；共享边（3 个点、中间是 `EDGE_PARALLEL`）返回 0（`common_subdivision.cpp:300-314`）。
- **[代码]** `getHalfedgePathAonB/BonA` 直接返回 `pointsAlong*` 中另一侧的 `pos*`（`common_subdivision.cpp:316-344`），**不滤除 `EDGE_PARALLEL`**；而整数坐标后端的 `trace*HalfedgeAlong*` 显式滤除了它（`integer_coordinates_intrinsic_triangulation.cpp:62, 77, 83, 97, 113, 119`）。直接遍历 `pointsAlong*` 的调用方若不过滤，会把 `t = 0.5` 这个伪中点当成真实交点。
- **[代码]** 共享边的恢复是启发式：在 A 边上看到两个**相邻的 `posB.type == Vertex`** 就插入 `EDGE_PARALLEL`，再在后者出发的半边中搜索指向前者的边（`signpost_intrinsic_triangulation.cpp:873-897`）。源码自称"90% sure"，以 `GC_SAFETY_ASSERT(found, ...)` 与 `pointsAlongB[eB].size() == 2` 兜底（897、904 行）。同一对顶点间存在多条边时不唯一。

### 2.4 面

**[代码]** `sourceFaceA` / `sourceFaceB` 把每个 CS 面映射到唯一包含它的 A 面 / B 面，`copyFromA/B` 因此逐面精确（`common_subdivision.ipp:28-46`）。CS 面由 `sliceFace` 按三条边上的**交点个数**切片（`common_subdivision.cpp:553-651`），只依赖组合信息，不依赖交点间的几何距离——所以组合结构对追踪误差稳健，几何位置完全由 `posA/posB` 决定。

### 2.5 哪些连续性得以保留

| 属性类型 | 通道 | 保留 | 丢失 |
| --- | --- | --- | --- |
| 顶点标量（PL） | Pointwise / L2 | 目标网格上的 C⁰；公共顶点值（Pointwise 精确） | A 边上的梯度跳变被摊成宽度约一条 B 边的斜坡；L2 把误差扩散到整个网格 |
| 面常量 | `copyFromA/B` | 逐 CS 面精确；跳变位置随 A 边映射，幅度不变 | 无 —— 这是唯一能精确表达"面内常数 + 跨边跳变"的通道 |
| 边常量 | 无 API | — | 需自行用 `pointsAlong*` + 几何重建 |
| 向量 / 方向 | Pointwise / L2（仅 double） | 逐分量数值 | **不做平行移动**（见 5.4） |

结构性限制[代码]：顶点 PL 函数在曲面上处处连续，所以"跨 A 边有跳变"的属性用顶点通道**无法表达**，只能被插值成斜坡；这不是数值误差，是函数空间的限制。只有面常量通道能保留跳变。

### 2.6 退化与多义情形：结果与保证边界

| 情形 | 源码行为 | 结果 | 边界 |
| --- | --- | --- | --- |
| A 顶点落在 B 边内部 | 抛 `"encountered vertex intersection in the middle of an edge"`（`common_subdivision.cpp:424-427`；`integer_coordinates_intrinsic_triangulation.cpp:363-369`） | 异常 | 硬失败，不静默 |
| B 顶点落在 A 边内部 | signpost 记 `EDGE_VERTEX` 并加入 `pointsAlongA`（791-798）；整数坐标记 `EDGE_VERTEX` 但不加入（213-216，改由法坐标追踪表示，238-349） | 正常 | 迁移经 `P_A/P_B` 时等价；遍历 `pointsAlongA` 时两后端给出不同点集 **[实测]** |
| signpost 把边上的点记成面点 | 见 2.2 | 数值正确，但 CS 含 2 个零面积面 | **无报错**；仅影响按面遍历的下游 **[实测]** |
| 两条边共线且部分重叠（T 形接合） | 切片与 `EDGE_PARALLEL` 恢复都假设"空 / 单交点 / 整条 B 边两端点为 A 顶点"三种之一（`signpost_intrinsic_triangulation.cpp:873-906`） | 组合上缺边界段 | **推断不受支持**；仓库无测试，未实测 |
| 退化面（边数 ≤ 2） | `std::cerr << "Error: degree-2 face in common refinement"; continue;`（`common_subdivision.cpp:455-464`） | 只打一行、跳过 | 若该点仍被其它面引用则 CS 有洞（M 行和 ≠ 总面积）；若无人引用则 `ManifoldSurfaceMesh` 抛 `"unreferenced vertex"`（`manifold_surface_mesh.cpp:145-149`） |
| 切片奇偶性不满足 | `GC_SAFETY_ASSERT((nij+njk+nki)%2==0, ...)`（`common_subdivision.cpp:622-625`） | 异常 | 硬失败 |
| 找不到 CS 面的父面 | `GC_SAFETY_ASSERT(sharedFace != Face(), ...)`（`common_subdivision.cpp:531`），注释承认策略"not necessarily an ideal strategy"（469-471） | 异常 | 默认构建硬失败；`NGC_SAFETY_CHECKS` 下 assert 被编译掉（`utilities.h:14-21`），随后拿到无效 `Face` |
| 交点 `orientation` | signpost 两处硬编码 `false`，注释 "TODO not properly populated"（`signpost_intrinsic_triangulation.cpp:840, 883`） | 恒为 false | 依赖"交点哪一侧"的逻辑在 signpost 后端不可用；整数坐标填了 `positiveOrientation`(288, 299) |
| 追踪被 `maxIters` 截断 | `options.maxIters = nFaces*10`（`signpost_intrinsic_triangulation.cpp:76`）；超限时把 `endPoint` 设在中途边上返回（`trace_geodesic.cpp:632-646`） | 追踪提前结束 | `recoverTraceTValues` 将累计长度**归一化到 [0,1]**（`intrinsic_triangulation.cpp:316-320`），整条边上所有交点被拉伸：组合结构仍对，**几何位置静默错误** |
| 追踪总长为 0 | `tVals[i] /= totalLen`（`intrinsic_triangulation.cpp:317-320`） | `0/0 = NaN` | **无保护**；NaN 的 `tEdge` 进入 `SurfacePoint`，而 `validate()` 不在迁移路径上调用（`surface_point.ipp:228-255`，全仓库仅 `flip_geodesics` 使用） |
| 点对点：零向量 | `traceGeodesic` 早退，`endPoint` 保持默认构造的**无效** `SurfacePoint`（`trace_geodesic.cpp:703-716`） | 无效点 | 无检测；`equivalentPointOn*` 先调用 `inSomeFace()` 任选一个相邻面（`surface_point.ipp:39-86`） |
| 整数坐标的点对点 | `equivalentPointOnIntrinsic` 直接 `throw "not implemented"`（`integer_coordinates_intrinsic_triangulation.cpp:166-170`） | 不可用 | 该后端**没有** A→B 点对应原语 |
| 最近的边吸附 | `TraceOptions` 中**并不存在** `allowEndOnEdge` / `allowEndOnEdgeEps` 字段（`trace_geodesic.h:22-27`），相关分支整段注释掉（`trace_geodesic.cpp:172-183, 396-407`）；`snapEndToEdgeIfClose` 仍是 TODO（`intrinsic_triangulation.cpp:274`；`signpost_intrinsic_triangulation.cpp:734`） | 未实现 | 退化情形本可用的兜底机制目前关闭 |

---

## 3. 问题二：三项目标是否兼容

### 3.1 形式化

设 A 的节点值 `f ∈ R^{V_A}`，B 的节点值 `g ∈ R^{V_B}`，`p = P_A f`，`q = P_B g`。"函数"指 CS 上的 PL 函数（`docs/.../function_transfer.md:26` 明确采用这一观点）。

- (T1) 节点复现：`g[k] = f[k]` 对所有共享节点 k。
- (T2) 往返恢复：A→B→A 后 `f` 逐位不变。
- (T3) 面积积分保持：`1ᵀM p = 1ᵀM q`。
- (T4) L2 最优：`‖p − q‖_M` 最小。

### 3.2 Pointwise ⇒ (T1)、(T2)[代码 + 实测]

**[代码]** `transferAtoB_Pointwise` 只在 `p.posB.type == Vertex` 时写值，值为 `p.posA.interpolate(f)`（`transfer_functions.cpp:63-69`）。对共享顶点 `posA = SurfacePoint(vA)`，`interpolate` 直接返回 `f[vA]`（`surface_point.ipp:203-206`）⇒ (T1)。反向同理（82-95）。由于构造时每个 B 顶点恰好压入一条记录（`signpost_intrinsic_triangulation.cpp:777-789`；`integer_coordinates_intrinsic_triangulation.cpp:199-212`），在"A 顶点在 B 中同索引存在"（§5.2 前提 2）下，A 的每个顶点在 CS 中恰有一个 `VERTEX_VERTEX` 点 ⇒ (T2)。

**[实测]** 插点与翻边两种情形下，A→B→A 都逐位返回原值；`sampleFromInput`→`restrictToInput` 亦然。

关于 (T3)：Pointwise 无任何保证。它只固定节点值，而两个三角划分对同一组节点值的插值积分一般不同（实测 `1/3 → 2/3`）。另外 Pointwise 保持值域（值为源节点值的凸组合），故不越界。

方向不对称[代码]：B→A→B **不是**恒等。`restrictToInput` 只抄写落在 A 顶点上的 B 顶点值（`intrinsic_triangulation.ipp:32-39`），插入顶点的值被丢弃、随后由 A 的 PL 函数重新插值。实测 `B = [0,1,2,3,4] → B' = [0,1,2,3,1]`。

### 3.3 L2 ⇒ (T3)、(T4)[代码 + 实测]

**[代码]** `transferAtoB_L2` 求解

```
(P_Bᵀ M P_B) g = P_Bᵀ M P_A f        (transfer_functions.cpp:72-80; 矩阵见 107-117)
```

即把 `p = P_A f` 做 **M-内积下的正交投影**到 `range(P_B)`：

- (T4)：残差与 `range(P_B)` M-正交 ⇒ 唯一最小解（`P_B` 列满秩时）。
- (T3)：`P_A`、`P_B` 每行都是重心权重且和为 1（`common_subdivision.cpp:110-136, 152-178` 的三个分支：`1`、`(1−t)+t`、`c0+c1+c2`），故 `P_A 1 = P_B 1 = 1`。取 `<p − q*, 1>_M = 0` 得

  ```
  1ᵀ M P_A f = 1ᵀ M P_B g     对一切 f 成立。
  ```

  等价算子形式 `1ᵀ M P_A = 1ᵀ M P_B T`。

**[实测]** 在两个后端（signpost / 整数坐标）上，插点情形下 `max|rA − rBᵀT| = 5.6e-17`（E2）；翻边正方形上 `1ᵀMP_A = [0.3333, 0.1667, 0.3333, 0.1667]`、`1ᵀMP_B = [0.1667, 0.3333, 0.1667, 0.3333]`（E1），与解析值逐位一致。

关键限定：该恒等式**只依赖"行和为 1"与"M 对称"**，不依赖对应位置的几何精度；它守恒的是"以 M 为测度的积分"，只有当 M 是真实面积矩阵时才是几何面积（见 5.3）。

### 3.4 冲突：可精确复算的反例

单位正方形，`A(0,0), B(1,0), C(1,1), D(0,1)`。网格 A 用对角线 AC（面 `ABC`、`ACD`），网格 B 用对角线 BD（面 `ABD`、`BCD`），即一次内蕴翻边；CS 为绕中心 `m` 的 4 个三角形（各 1/4）。

取 A 节点值 `f = [A:0, B:1, C:0, D:1]`（跨 AC 有梯度跳变）。**[实测]** 结果（附有理数校验，见附录 A）：

| 量 | 值 |
| --- | --- |
| `∫_A f` | `1/3` |
| Pointwise A→B 节点值 | `[0, 1, 0, 1]`（满足 T1） |
| Pointwise A→B 后 `∫_B` | **`2/3`**（T3 失败，+1/3） |
| L2 A→B 节点值 | `[0, 0.5, 0, 0.5]`（T1 失败） |
| L2 A→B 后 `∫_B` | **`1/3`**（T3 成立） |
| L2 的 `‖p−q*‖_M²` | `1/24` |
| Pointwise 的 `‖p−q‖_M²` | `1/6`（> 1/24，符合 T4） |
| L2 往返 A→B→A | `[0.25, 0.5, 0.25, 0.5]`，`‖Δ‖_M² = 1/32`（T2 失败） |
| Pointwise 往返 A→B→A | `[0, 1, 0, 1]`，逐位相等（T2 成立） |
| L2 把 B 上单位脉冲映回 A | `[0.25, 1, 0.25, ~0]`（无越界，见下） |
| L2 把 A 上单位脉冲（节点 1）映到 B | `[0, 0.75, 0, **−0.25**]`（越界） |

冲突的本质：两侧 PL 空间都是 CS 上 5 维空间的 4 维子空间，交集为 3 维。`f` 不在交集中，于是

- (T1)+(T2) 强制 `q = P_B f`（把同一组节点值放进 B 的空间），积分由 B 的三角划分决定，一般 ≠ `∫p`；
- (T3)+(T4) 强制 `q` 为正交投影，公共节点值被改动。

只有当 `p ∈ range(P_B)`（`f` 在两个三角划分的叠加单元上分片线性）时才同时成立。对平面四边形，叠加的 4 个单元上分片线性 ⇒ 该四边形上**仿射**；所以一次翻边时交集恰为四边形上的仿射函数（本例 3 维）。

**[实测]** 验证：`f = y`（`[0,1,1,0]`）时 Pointwise 与 L2 都精确复原、积分 `1/2` 一致、往返精确；`f = const 3` 同理。

补充性质[代码+实测]：A→B→A 的 L2 复合**不是幂等投影**（两个不同子空间的投影复合）。实测反复往返单调收敛到交集上的投影（常数 `1/3`），且每步积分恒为 `1/3`：

```
start        [0, 1, 0, 1]        mass 0.333333
round 1   [0.25, 0.5, 0.25, 0.5] mass 0.333333
round 2   [0.3125, 0.375, ...]   mass 0.333333
round 5   [0.333008, 0.333984, ...] mass 0.333333
```

所以"一次往返不精确"是定义使然，反复往返不会回到原值。

### 3.5 现有方案选择了哪些保证

| 目标 | `Pointwise` | `L2` | 依据 |
| --- | --- | --- | --- |
| 共享节点值复现 (T1) | ✅ 精确 | ❌ | `transfer_functions.cpp:63-69` / 72-80；实测 |
| A→B→A 恒等 (T2) | ✅ 逐位 | ❌ | 3.2 / 3.4 |
| 保积分 (T3) | ❌ | ✅（以所用 M 为测度） | 3.3 正交性；实测逐列验证 |
| 最小 L2 误差 (T4) | ❌ | ✅ 唯一 | 正规方程；测试 `intrinsic_triangulation_test.cpp:395` 断言 `L2Err ≤ PointwiseErr` |
| 值域 / 无过冲 | ✅ | ❌（实测 `−0.25`） | 凸组合 / 实测 |
| 与 B 的划分无关（节点值） | ✅ | ❌ | Pointwise 只依赖 `posA`；L2 依赖整个 `P_B` |
| 成本 | O(V_B)，无需 CS | 需 CS + 预分解 + 一次求解 | `transfer_functions.cpp:17-26` |

`AttributeTransfer` 提供 L2 的预分解（`AtoB_L2_Solver` / `BtoA_L2_Solver`，`transfer_functions.h:33-34`），重复迁移多组属性时可摊薄开销；代价是把"同一 M、同一 CS"的约束固化在对象生命周期里（见 5.2 前提 4）。

---

## 4. 实测证据（原始输出见第 6 节）

| 编号 | 探针 | 关键数值 |
| --- | --- | --- |
| E1 | 翻边正方形上的 Pointwise / L2 / 往返 | 见 3.4 表 |
| E2 | 逐列积分恒等式（两个后端，插点情形） | `max\|rA − rBᵀT\| = 5.6e-17`（两者） |
| E3 | 反复 L2 往返 | 收敛到常数 `1/3`，积分恒 `1/3` |
| E4 | signpost vs 整数坐标：共享内边插点（各一次独立运行） | `FACE_VERTEX`(posA 面内) + 额外 `EDGE_TRANSVERSE` → 6V6F，2 个零面积面 vs `EDGE_VERTEX`(t=0.5) → 5V4F，无退化面 |
| E5 | 针状三角形 h=1e-3，mollification | 边长加 `~4.7e-6`；`AttributeTransfer(intTri)` 面积 `9.13e-4`，真实 `5.0e-4` |
| E6 | 同一针状体，构造后 `faceAreas` | 报 `5.0000000001e-4`（mollification 前），`purgeQuantities()` 后变 `9.13e-4` |
| E7 | `AttributeTransfer` 跨变异持有 | L2 路径**段错误**（exit 139，访问已释放的 CS）；Pointwise 路径静默返回旧规模的 `P_A`（4 行 vs 当前 CS 5 顶点） |
| E8 | `restrictToInput` / B→A→B | 前者精确；后者插入顶点值丢失 `[0,1,2,3,4] → [0,1,2,3,1]` |
| E9 | 弯曲"帐篷"网格上 M 的两种构造 | `‖M_A − M_B‖ = 0` |
| E10 | 平坦正方形上 M 的两种构造 | `‖M_A − M_B‖ ≈ 1.8e-16` |
| E11 | 仿射函数 `f = x` 在两个后端上的 L2 迁移 | 最大误差 `1.7e-16`；两后端数值一致 |
| E12 | 零面积面的稳定性 | `t = .25/.3/.5/.7` 均为 6V6F + 2 个零面积面；闭合四面体（flip+split）为 6V8F + 0 个 |
| E13 | 针状体上两种 M 构造（整数坐标，moll=1e-5） | `M(positionsA)` 总面积 `5.000e-4`（=真值），`M(lengthsB)` 为 `9.129e-4`，`‖M_A−M_B‖ = 1.46e-4` |

---

## 5. 问题三：选择条件与结果判据

### 5.1 何时选哪条通道

| 需求 | 选择 | 理由 / 注意 |
| --- | --- | --- |
| 有界/非负属性（浓度、指示函数、不透明度），要求节点值逐位一致 | Pointwise | 凸组合保界；A→B→A 精确 |
| 关心积分/总量/平均，允许节点值改变 | L2 | 正交投影；须确认 M 是真实面积矩阵（5.3） |
| 既要节点一致又要积分一致 | 仅当函数属于两侧空间交集（叠加单元上分片线性，如纯插点或仿射） | 否则无解；用 `‖p−q‖_M` 与 `1ᵀM(p−q)` 两个判据同时检验，后者 ≠ 0 即已选择 L2 语义 |
| 面内整数标签/材质/面法向 | `copyFromA/B` | 逐面精确，保留跳变位置与幅度 |
| 边常量 | 无现成 API | 需从 `pointsAlong*` + 几何重建，并滤除 `EDGE_PARALLEL` |
| 向量/方向场 | 两者都不安全 | 需显式平行移动；`AttributeTransfer` 仅 `double` |
| 只需可视化 / 在 CS 上运算 | `interpolateAcrossA/B` | 升到 CS，不做折损；注意 signpost 的零面积单元（2.2） |

### 5.2 结果判据（可直接写成断言）

1. **节点恒等**：`restrictToInput(sampleFromInput(f))[v] == f[v]` 对 A 的每个顶点逐位成立。B→A→B 只应在公共顶点上成立。
2. **正规方程残差**：`‖(P_BᵀMP_B)g − P_BᵀMP_A f‖ ≤ 1e-6`（同 `intrinsic_triangulation_test.cpp:397-400`），或用 `constructAtoBMatrices()`。
3. **投影残差下界**：`‖p − q*‖_M² == ‖p‖_M² − ‖q*‖_M²`；不等说明求解异常或 M 不对称。
4. **积分差**：`1ᵀMP_A f − 1ᵀMP_B g`。L2 应为 0（容差内）；Pointwise 一般非 0，可量化"选了 Pointwise 语义"。
5. **包络**：`min(g) ≥ min(f)`、`max(g) ≤ max(f)`。违反则来自 L2 且属性受有界约束，应改 Pointwise 或投影回可行域。
6. **度量一致性**：`cs.interpolateEdgeLengthsA(geomA.edgeLengths)` 与 `...B(intTri.edgeLengths)` 之差（`intrinsic_triangulation_test.cpp:432-434` 的做法）；`vertexGalerkinMassMatrixFromPositionsA` 与 `FromLengthsB` 的行和之差（仓库**没有**这条测试）。
7. **CS 完整性**：`cs.mesh->nFaces()` 与 `cs.elementCounts()` 一致（`intrinsic_triangulation_test.cpp:309-315`）；stderr 无 `"degree-2 face"`；**无零面积面**（signpost 会引入，见 2.2——仓库没有这条检查）。
8. **追踪健全性**：每条 B 边追踪路径的首/末点应对应两个端点；`recoverTraceTValues` 严格单调且 `totalLen > 0`。
9. **点对应**（仅 signpost）：`‖pos(equivalentPointOnInput(equivalentPointOnIntrinsic(p))) − pos(p)‖ < 1e-5`（`intrinsic_triangulation_test.cpp:115-127`）。

### 5.3 哪些保证由源码支持，哪些依赖调用方前提

**源码支持**

- **[代码]** `P_A`、`P_B` 行和为 1 ⇒ 常量精确保持；这是 L2 积分恒等式的支点（`common_subdivision.cpp:110-178`）。
- **[代码]** Pointwise 在公共顶点逐位复现，A→B→A 逐位恒等；B→A→B 有损（3.2）。
- **[代码]** L2 是 M-正交投影 ⇒ 唯一最小 L2 误差 + 对任意数据保持 `1ᵀM` 积分（3.3）。
- **[代码]** `L2Err ≤ PointwiseErr` 由定义成立，非经验结论。
- **[代码]** 纯欧氏翻边保持度量：新边长由 `layoutDiamond` 在平面菱形上量得（`signpost_intrinsic_triangulation.cpp:289-294, 330-343`）；插点同理。度量保真时两种 M 构造一致（**E9/E10 实测** `‖M_A−M_B‖ ≤ 2e-16`）。
- **[代码]** 边界边不可翻（`isFixed` 对 `isBoundary` 返回 true，`intrinsic_triangulation.ipp:8-12`），故 B 与 A 的组合边界一致，边界属性可按一维插值精确迁移。
- **[代码]** 任何变异都清空 CS（`triangulationChanged()` → `commonSubdivision.reset()`，`intrinsic_triangulation.cpp:712`；调用点 `signpost_intrinsic_triangulation.cpp:319, 370, 395, 480, 532, 576`；`integer_coordinates_intrinsic_triangulation.cpp:433, 503, 1407, 1561, 1681, 1705`）⇒ "变异后必须重建 CS"是设计契约。

**调用方必须提供的前提**

1. **[代码]** A 已 `compress()` 且为三角网格（构造器检查并抛异常，`intrinsic_triangulation.cpp:19-31`）。
2. **[代码]** 每个 A 顶点在 B 中以相同索引存在（`docs/.../common_subdivision.md:11`；用 `getIndex()` 对应）。这不是运行时检查项：整数坐标后端若对应被破坏，只会得到错误位置。原始顶点永不被删、始终占据 `[0, n_A)`，由 `compress()` 保序（`surface_mesh.cpp:1900-1930`）加上 `removeInsertedVertex` 拒绝删除原始顶点（`signpost_intrinsic_triangulation.cpp:545`；`integer_coordinates_intrinsic_triangulation.cpp:1585`）共同维持，属实现细节而非文档承诺。
3. **[代码]** 在 `IntrinsicTriangulation` 生命周期内不修改 `inputMesh` / `inputGeom`（`docs/.../basics.md:105, 110`）。
4. **[代码+实测 E7]** 持有 `AttributeTransfer` 期间不要变异内蕴三角划分。`AttributeTransfer::cs` 是对 `IntrinsicTriangulation` 所拥有 CS 的**引用**（`transfer_functions.h:27`），构造器把 `P_A`、`P_B`、`M_CS_Galerkin` 拷出（`transfer_functions.cpp:12-25`），而变异会 `triangulationChanged()` → `commonSubdivision.reset()`（`intrinsic_triangulation.cpp:712`）删除被引用的对象。**实测**：在 `transferAtoB(..., L2)` 之后翻一条边，再对本对象调用 `transferAtoB(..., L2)`，进程**段错误（exit 139）**——L2 路径上的 `AtoB_L2_Solver` 已按旧 CS 分解，随后访问已释放的 `cs`/`mesh`。Pointwise 路径不触碰 `cs`，实测仍返回旧 CS 规模的结果（`P_A` 行数 4 vs 当前 CS 顶点 5），即**静默错误**而非崩溃。文档只提醒"如需保留请复制 CS"（`docs/.../basics.md:147-153`），构造函数本身不拷贝、也不注册任何失效回调。
5. **[代码+实测 E5/E6]** 用 L2 追求"面积保持"时，M 必须是真实面积矩阵：
   - `AttributeTransfer(CommonSubdivision&, VertexPositionGeometry&)` 用 A 的顶点位置构造 M（`transfer_functions.cpp:14`；面积由叉积计算，`embedded_geometry_interface.cpp:275-292`）；
   - `AttributeTransfer(IntrinsicTriangulation&)` 用 B 的边长构造 M（`transfer_functions.cpp:25`；Heron，`edge_length_geometry.ipp:9-26`）；
   - 度量保真且 CS 完整时二者相同（实测 E9/E10 `‖M_A−M_B‖ ≤ 2e-16`）。
   - **但整数坐标后端默认 mollification**：`mollifyIntrinsic(*intrinsicMesh, edgeLengths, 1e-5)`（`integer_coordinates_intrinsic_triangulation.cpp:33-41`；`intrinsic_mollification.cpp:20-38` 给所有边加同一个正常数，只要存在 `lC > lA + lB` 的退化面）。实测 E5/E13：h=1e-3 的针状三角形上，`M(positionsA)` 总面积 `5.000e-4`（=真值），`M(lengthsB)` 为 `9.129e-4`，即 **1.83 倍**，`‖M_A−M_B‖ = 1.46e-4`。此时 L2 的"面积守恒"守恒的是错误的面积；用 A 构造的迁移与用 B 构造的迁移不再互为伴随。
   - 且 mollification **发生在 `faceAreas` 已被缓存之后**：基类构造器 `requireHalfedgeVectorsInFace()`（`intrinsic_triangulation.cpp:50`）→ `faceAreasQ.ensureHave()`（`intrinsic_geometry_interface.cpp:292`）→ 计算 `faceAreas`；随后 mollify 就地改 `edgeLengths`（第 39-41 行）而无失效。实测 E6：构造后 `faceAreas` 仍是 mollification 前的值，`purgeQuantities()` 后才更新。于是同一个 `IntrinsicTriangulation` 上 `edgeLengths` 与 `faceAreas` 不一致。
6. **[代码]** 属性是标量，或与切标架无关（见 5.4）。
7. **[代码]** 不依赖 signpost 后端的 `orientation`（恒为 false）。
8. **[代码]** 自行遍历 `pointsAlong*` 时需滤除 `EDGE_PARALLEL`（`t = 0.5` 伪点）。

### 5.4 向量/方向场的缺口

**[代码]** `AttributeTransfer` 只接受 `VertexData<double>`（`transfer_functions.h:40-52`），无向量通道。`sampleFromInput/restrictToInput` 虽模板化，但 `SurfacePoint::interpolate` 对标量/向量一视同仁做逐分量线性插值（`surface_point.ipp:200-226`），不做平行移动。而顶点切基由 `signpostAngle` / `halfedgeVectorsInVertex` 给出（`signpost_intrinsic_triangulation.cpp:589-625`；`intrinsic_geometry_interface.cpp:367-392`），其 `v.halfedge()` 在翻边后可能改变，因此同一顶点在 A、B 中的基可以相差一个旋转。文档示例却直接用 `restrictToInput` 搬运方向场（`docs/.../basics.md:59`；`docs/.../algorithms/vector_heat_method.md:126`），同时承认切空间对齐"may or may not be a builtin method"（该文件 129 行）。

---

## 6. 尚缺的证据

以下各项在仓库中没有测试或实现覆盖（已用 `grep` 确认，命名为便于补测）：

1. **`sampleFromInput` / `restrictToInput` 零测试**：`grep` 在整个 `test/` 中命中 0 次，覆盖全部经由 `AttributeTransfer`。
2. **无"只翻边、不插点"的真值测试**。`FunctionTransfer`（`intrinsic_triangulation_test.cpp:373-402`）用整数坐标 + `delaunayRefine()`（含 mollification），输入是 `Vector<double>::Random`，只断言 `L2Err ≤ PointwiseErr` 与残差 ≤ 1e-6。前者由最优性自动成立，后者只验证"求解器解了自己的方程"——**M 错了也会通过**。
3. **无积分守恒/不守恒断言、无有界性断言、无"L2 改变公共节点值"的断言**。3.4 的反例不在仓库中。
4. **无开边界网格的迁移测试**，尽管边界保持可由源码证明。
5. **无退化 CS 测试**：`degree-2` 跳过（仅 stderr）、部分共线重叠、重边/多平行边、`VERTEX_VERTEX` 落在边内。`CommonSubdivisionCompareIntegerSignpost`（328-368）只比较元素计数。
6. **无零面积 CS 面的检查**——signpost 后端在开边界网格的内边上插点时稳定产生 2 个（实测 E4/E12），而库只用 `"degree-2 face"` 覆盖"边数 ≤ 2"，不覆盖"面积 = 0 的三角形"。
7. **无 `M_CS_Galerkin` 两种构造器 / 两个后端的一致性测试**。`CommonSubdivisionGeometry`（404-436）只比较 CS 的**边长**，不比较质量矩阵或积分。
8. **无 mollification 对迁移影响的测试**，也无"mollification 后 `faceAreas` 缓存失效"的检查（实测 E5/E6 暴露的实际缺陷）。
9. **无追踪失败率、`maxIters` 截断、`trimTraceResult` 失败的统计或断言**；`recoverTraceTValues` 的 `totalLen == 0` 无保护也无测试。
10. **无"`AttributeTransfer` 跨变异"的测试或保护**：实测 L2 路径段错误、Pointwise 路径静默错误（E7）。构造函数与 `transferAtoB_L2` 都没有任何 CS 有效性校验。
11. **`SquareSolver` 在 SuiteSparse 路径下不检查分解状态**：`umfFactor` 返回 `void`，`umfpack_*` 的返回值被忽略（`square_solvers.cpp:53-59, 129-139`）；Eigen 路径会抛异常（144-147）。CS 退化导致的奇异系统在 SuiteSparse 构建下可能给出静默错误解。本机构建无 SuiteSparse，故**未实测**该分支。
12. **边界情形：两条边共线部分重叠**——本文件的"不受支持"是**代码分支推断**，未构造实例实测；需要单独验证。
13. **`equivalentPointOn*` 的歧义/退化无测试**：零向量早退返回无效点（`trace_geodesic.cpp:703-716`）、cut locus 多解、边界顶点角度归一化（`signpost_intrinsic_triangulation.cpp:195-201`）都只有"经验上可行"的注释。
14. **mollification 的对象一致性问题**（E5/E6/E13）需要单独修复验证：`faceAreas` 在 mollification 之后未失效（实测需 `purgeQuantities()` 才更新），且 `mollifyEPS` 只改 `intrinsicMesh` 的边长、不改 `inputMesh`/`geomA`。实测 E13：同一针状体上 `M(positionsA)` 总面积 `5.000e-4`（=真值），`M(lengthsB)` 总面积 `9.129e-4`，`‖M_A − M_B‖ = 1.46e-4`。仓库无此测试。

---

## 附录 A：3.4 反例的精确复算（有理数，无浮点）

```python
from fractions import Fraction as F
# 单位正方形 0=A(0,0) 1=B(1,0) 2=C(1,1) 3=D(0,1)
#   A: 对角线 AC -> (A,B,C),(A,C,D)；B: 对角线 BD -> (A,B,D),(B,C,D)
#   CS: 中心 m 周围 4 个三角形，面积各 1/4
faces=[[0,1,4],[1,2,4],[2,3,4],[3,0,4]]; areas=[F(1,4)]*4; n=5
M=[[F(0)]*n for _ in range(n)]
for f,a in zip(faces,areas):
    for i in range(3):
        M[f[i]][f[i]]+=a/6; M[f[i]][f[(i+1)%3]]+=a/12; M[f[i]][f[(i+2)%3]]+=a/12
PA=[[F(1),0,0,0],[0,F(1),0,0],[0,0,F(1),0],[0,0,0,F(1)],[F(1,2),0,F(1,2),0]]
PB=[[F(1),0,0,0],[0,F(1),0,0],[0,0,F(1),0],[0,0,0,F(1)],[0,F(1,2),0,F(1,2)]]
T   =lambda P,v:[sum(P[r][c]*v[c] for c in range(len(v))) for r in range(len(P))]
mv  =lambda A,v:[sum(A[i][j]*v[j] for j in range(len(v))) for i in range(len(v))]
dot =lambda a,b: sum(x*y for x,y in zip(a,b))
w   =mv(M,[F(1)]*n); integ=lambda v: dot(w,v)
wPA =[sum(w[i]*PA[i][j] for i in range(n)) for j in range(4)]
wPB =[sum(w[i]*PB[i][j] for i in range(n)) for j in range(4)]
def gauss(A,b):
    N=len(A); Tm=[row[:]+[b[i]] for i,row in enumerate(A)]
    for c in range(N):
        p=next(r for r in range(c,N) if Tm[r][c]!=0); Tm[c],Tm[p]=Tm[p],Tm[c]
        pv=Tm[c][c]; Tm[c]=[x/pv for x in Tm[c]]
        for r in range(N):
            if r!=c and Tm[r][c]!=0:
                fa=Tm[r][c]; Tm[r]=[x-fa*y for x,y in zip(Tm[r],Tm[c])]
    return [Tm[i][N] for i in range(N)]
def l2op(Ps,Pt,ns):
    oh=lambda j:[F(1) if k==j else F(0) for k in range(ns)]
    it=[T(Pt,oh(j)) for j in range(ns)]; isrc=[T(Ps,oh(j)) for j in range(ns)]
    A=[[dot(it[j],mv(M,it[k])) for k in range(ns)] for j in range(ns)]
    cols=[gauss(A,[dot(it[i],mv(M,isrc[j])) for i in range(ns)]) for j in range(ns)]
    return [[cols[j][i] for j in range(ns)] for i in range(ns)]
TAB=l2op(PA,PB,4); TBA=l2op(PB,PA,4)
ap=lambda Tm,v:[sum(Tm[i][j]*v[j] for j in range(len(v))) for i in range(len(Tm))]
f=[F(0),F(1),F(0),F(1)]
print("∫_A f =", integ(T(PA,f)), " ∫ pointwise =", integ(T(PB,f)))
g=ap(TAB,f); print("L2 g =",g," ∫ =",integ(T(PB,g)))
h=ap(TBA,g); d=[T(PA,h)[i]-T(PA,f)[i] for i in range(n)]
print("L2 往返 =",h," ‖Δ‖_M² =",dot(d,mv(M,d)))
print("L2 脉冲越界 =", ap(TAB,[F(0),F(1),F(0),F(0)]))
print("wPA =",wPA," wPB =",wPB)
print("行恒等 =", wPA==[sum(wPB[i]*TAB[i][j] for i in range(4)) for j in range(4)])
```

输出：

```
∫_A f = 1/3  ∫ pointwise = 2/3
L2 g = [0, 1/2, 0, 1/2]  ∫ = 1/3
L2 往返 = [1/4, 1/2, 1/4, 1/2]  ‖Δ‖_M² = 1/32
L2 脉冲越界 = [0, 3/4, 0, -1/4]
wPA = [1/3, 1/6, 1/3, 1/6]  wPB = [1/6, 1/3, 1/6, 1/3]
行恒等 = True
```

**[实测]** 用本仓库真实 API 重跑同一构型（`SignpostIntrinsicTriangulation` + `flipEdgeIfPossible` + `AttributeTransfer`，无 SuiteSparse）：

```
flipped = 1  B verts=4 faces=2 edges=5
int_CS(P_A f)     = 0.333333
pointwise B       = 0 1 0 1  int=0.6666666667  err2=0.1666666667
L2        B       = 0 0.5 0 0.5  int=0.3333333333  err2=0.04166666667
PW  A->B->A       = 0 1 0 1
L2  A->B->A       = 0.25 0.5 0.25 0.5
1^T M P_A = 0.333333 0.166667 0.333333 0.166667
1^T M P_B = 0.166667 0.333333 0.166667 0.333333
```

与解析解逐位一致。

## 附录 B：复现方式

```bash
# 依赖：Eigen 3.4.0（解压到 /tmp/gc_deps/eigen-3.4.0）
cmake <repo> -DGC_EIGEN_LOCATION=/tmp/gc_deps/eigen-3.4.0 -DCMAKE_BUILD_TYPE=Release
cmake --build . -j8
# 手写 main() 链接 build/src/libgeometry-central.a 即可复现第 6 节各探针
```

## 附录 C：关键源码位置索引

| 主题 | 位置 |
| --- | --- |
| CS 数据结构与交点类型 | `include/geometrycentral/surface/common_subdivision.h:18-44, 100-153` |
| CS 构造（signpost） | `src/surface/signpost_intrinsic_triangulation.cpp:762-909` |
| CS 构造（整数坐标） | `src/surface/integer_coordinates_intrinsic_triangulation.cpp:187-372` |
| 切片 / 退化面跳过 / 奇偶性断言 | `src/surface/common_subdivision.cpp:455-464, 538-551, 553-651` |
| 插值矩阵 | `src/surface/common_subdivision.cpp:100-182` |
| CS 边长与质量矩阵 | `src/surface/common_subdivision.cpp:184-256`；`src/surface/intrinsic_geometry_interface.cpp:470-501` |
| 逐点迁移 / L2 迁移 | `src/surface/transfer_functions.cpp:57-117` |
| 逐点采样 / 限制 | `include/geometrycentral/surface/intrinsic_triangulation.ipp:22-40` |
| 追踪 t 值恢复 | `src/surface/intrinsic_triangulation.cpp:293-322` |
| 点对点对应（signpost） | `src/surface/signpost_intrinsic_triangulation.cpp:149-265` |
| 点对点对应（整数坐标，未实现） | `src/surface/integer_coordinates_intrinsic_triangulation.cpp:166-183` |
| 追踪器与早退 | `src/surface/trace_geodesic.cpp:632-646, 687-745` |
| mollification | `src/surface/intrinsic_mollification.cpp:20-38`；`src/surface/integer_coordinates_intrinsic_triangulation.cpp:33-41` |
| 线性求解器 | `src/numerical/square_solvers.cpp:117-187` |
| 迁移语义文档 | `docs/docs/surface/intrinsic_triangulations/function_transfer.md:1-52, 98-106` |
| CS 前置条件与几何一致性 | `docs/docs/surface/intrinsic_triangulations/common_subdivision.md:11, 103-113, 166-170` |
| 迁移相关测试 | `test/src/intrinsic_triangulation_test.cpp:105-153, 275-324, 328-368, 373-402, 404-436` |
