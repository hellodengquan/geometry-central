# 三角网格测地线追踪：穿边、过顶点、命中边界时的终点语义

分析对象：`src/surface/trace_geodesic.cpp`（`traceGeodesic` 的全部实现）、`include/geometrycentral/surface/trace_geodesic.h`、`include/geometrycentral/surface/surface_point.h`、`include/geometrycentral/surface/barycentric_coordinate_helpers.ipp`、`src/surface/intrinsic_geometry_interface.cpp`（切空间构造），以及调用方 `src/surface/signpost_intrinsic_triangulation.cpp`、`src/surface/intrinsic_triangulation.cpp`、`src/surface/integer_coordinates_intrinsic_triangulation.cpp`。

下文所有结论均标注源码位置。凡是通过实际编译运行验证的，标为「实测」并给出观测值（用仓库源码 + Eigen 3.4 构建，网格为平面/折面/锥点小网格）。无法从代码推出保证的地方明确标为「无保证」。

---

## 0. 结论速览

| 问题 | 实际答案 |
|---|---|
| 直行是曲面内蕴的还是三维的？ | **纯内蕴**。只依赖边长度，完全不读三维坐标（`trace_geodesic.cpp:689-691` 只 require `vertexAngleSums` / `halfedgeVectorsInVertex` / `halfedgeVectorsInFace`）。三维折线长度 == 报告长度，但三维路径在每条棱处有折角，不是三维直线。 |
| 正常耗尽长度时终点是什么？ | `SurfacePointType::Face`，`hitBoundary=false`，`length == norm(traceVec)`。 |
| 命中边界/barrier 时终点是什么？ | `SurfacePointType::Edge`，`hitBoundary=true`，`length < norm(traceVec)`（实际走过长度）。所在边可能是真正的边界边，也可能是内部 barrier 边。 |
| 精确命中内部边呢？ | **通常被算作 Face**，不是 Edge（`isInsideTriangle` 用 `>= 0`，`barycentric_coordinate_helpers.ipp:70`）。实测 16/16 为 Face。但本质上是浮点结果决定，与下面边界边的情形同理。 |
| 精确命中边界边呢？ | **数值上不确定**：可能 Face+`hitBoundary=false`，也可能 Edge+`hitBoundary=true`，取决于浮点结果是恰好 0 还是略负。 |
| 过顶点时终点是什么？ | 若恰好耗尽在顶点：`Face` 类型但 `faceCoords=(1,0,0)`，几何上就是顶点，**没有任何「这是顶点」的标记**。若越过顶点继续：会产生一串 `t=0`/`t=1` 的 Edge 点（3D 上全是同一个点），然后继续。**永不产生 Vertex 类型的中间点**。 |
| `maxIters` 管什么？ | 管「跨边 + 走完下一个面」的步数，不是面数上限也不是长度。`maxIters=0` 表示停在**第一个交点**。截断时 `hitBoundary` 保持 `false`，且终点必是内部非 barrier 边。 |
| 长度怎么记账？ | `length` 先置为 `norm(traceVec)`，终止时减去「剩余长度」。正常结束减 0；边/barrier/截断减 `lenRemaining`。 |
| 已知 bug | 从边界边上、朝外追踪会**立即终止但不清算长度**，`TraceSubResult::traceVectorInHalfedgeLen` 未初始化（`trace_geodesic.cpp:493-501`），`length` 读到栈上垃圾值（实测污染栈后稳定为 `-7.7e307`）。 |
| `endingDir` 的坐标系？ | 文档说是「该点的切空间」（即 `edge.halfedge()` 方向为 x 轴），但**边上的终止实际用的是「被跨过的那条半边 `crossHe`」的坐标系**。网格边界边的停止恰好与规范坐标系一致；但**内部 barrier 边的停止有约一半情况相差 180°**（实测 8/16 例翻转）。 |

---

## 1. 数据流与两种几何

### 1.1 入口与依赖量

两个重载（`trace_geodesic.cpp:687`、`:751`）都只 require 三个内蕴量：

```cpp
geom.requireVertexAngleSums();
geom.requireHalfedgeVectorsInVertex();
geom.requireHalfedgeVectorsInFace();
```

三者全部由 `edgeLengths` 推出（`intrinsic_geometry_interface.cpp:135-145`、`:367-392`、`:291-339`），不涉及顶点坐标。**实测**：同一网格分别用 `VertexPositionGeometry` 和由同一组边长构造的 `EdgeLengthGeometry` 追踪，沿 3 个方向得到的 `endPoint`、`pathPoints`、`length` 逐位相同（差异 0）。这就是「内蕴」的结构性证据。

### 1.2 三种起点

`traceGeodesic(SurfacePoint, Vector2)` 按 `startP.type` 分派（`:723-738`）：

- `Vertex` → `traceGeodesic_fromVertex`（`:523`）
- `Edge` → `traceGeodesic_fromEdge`（`:473`）
- `Face` → `traceGeodesic_fromFace`（`:510`）

另有 `traceGeodesic(Face, Vector3 startBary, Vector3 traceBaryVec)`（`:751`），先把 `startBary` 投影进三角形（`:788`），把 `traceBaryVec` 平移成和为 0 的位移（`:789`），再转成笛卡尔长度（`:793-794`）。

### 1.2.1 同一几何位置、不同类型 → 完全不同的切空间

`SurfacePointType` 不只是「描述这个点长什么样」，它**决定用哪个切空间解释输入方向**。最尖锐的例子是「几何上位于内部顶点、但类型是 Face」的点（`faceCoords=(1,0,0)`）：

- 类型 `Vertex`：`traceGeodesic_fromVertex` 走顶点扇区逻辑，x 轴 = `vertex.halfedge()` 方向，角度按 `angleSum/2π` 缩放（§4.1）。
- 类型 `Face`：`traceGeodesic_fromFace` 直接在**该面的二维坐标**里沿输入方向走，x 轴 = **`face.halfedge()` 的方向**，无角度缩放（`:510-519`）。

这两个切空间彼此相差一个旋转（以及可能的缩放）。**实测**（平面网格，同一个几何点 v4，同样输入 `0.2 * fromAngle(θ)`，θ 取 12 个值）：两者的三维位移**没有一次相同**（差异远超 1e-9）；例如 θ=0 时 Face 起点给出 `(0.2000, 0.0000)`，Vertex 起点给出 `(0.1414, 0.1414)`。

> 实践含义：**构造起点时必须让类型与几何位置一致。**把顶点写成 `SurfacePoint(face, Vector3{1,0,0})` 是合法的、不会报错，但会用错切空间。同理，`SurfacePoint(edge, t)` 与 `SurfacePoint(face, bary)` 即使几何重合，方向解释也不同。

### 1.3 三种行走子过程

| 子过程 | 位置 | 起点 | 允许跨越的边 |
|---|---|---|---|
| `traceInFaceBarycentric` | `:132` | 一般面内点 | 由 `edgeIsHittable` 指定 |
| `traceInFaceTowardsEdge` | `:328` | 顶点（指向某个边） | 只允许对边 `hittable[(iHe+1)%3]`（`:379-380`） |
| `traceInFaceFromEdge` | `:391` | 边上一点，进入相邻面 | 三条边都允许，唯独进入时所在边不允许（`:464-465`） |

主循环 `traceGeodesic_iterative`（`:626`）反复调用 `traceInFaceFromEdge`，每步「跨一条边 + 走完下一个面」。

---

## 2. 单面内的终止 / 跨越判据

`traceInFaceBarycentric` 是唯一决定单面内行为的地方。设面内位移 `vecBary`，终点 `endPoint = startPoint + vecBary`（`:162`）。

### 2.1 终止在面内（长度耗尽）

```cpp
if (isInsideTriangle(endPoint)) {            // :199
  result.terminated = true;
  result.endPoint = SurfacePoint(face, endPoint);   // Face 类型
  result.incomingDirToPoint = vecCartesianDir;      // 单位方向，面切空间
  result.traceVectorInHalfedgeLen = 0;              // :205  剩余长度清零
  return result;
}
```

`isInsideTriangle` 要求三个重心坐标都 `>= 0`（`barycentric_coordinate_helpers.ipp:70`）。**关键在于它是闭集判据**：重心坐标恰好为 0（点落在边或顶点上）也算「在面内」。

由此得到三条重要推论：

1. **精确落在内部边上 → 通常报 `Face` 类型，不是 `Edge`。** 该点的 `faceCoords` 有一个分量为 0。判据完全取决于浮点结果恰好是 `0.0` 还是略负，代码里没有任何「吸附到边」的逻辑。
   **实测**：对 8 个面的 16 条内部边「从重心精确瞄准边中点、长度恰好够到」，结果 16/16 全是 `Face`，`hitBoundary=false`。换成别的起点/方向，同样的瞄准方式也可能落到另一分支。
2. **精确落在内部顶点上 → 报 `Face` 类型，`faceCoords` 形如 `(1,0,0)`。** 几何上是顶点，但类型上完全看不出。
   **实测**：平面网格面 0 内从 `(0.25,0.25,0.5)` 精确瞄准顶点 v4，长度恰好等于距离 → `endType=Face`、`faceCoords=(1.000,0.000,0.000)`、`hitBoundary=false`。
3. **精确落在边界边上 → 结果数值上不确定。** 若重心坐标恰好算出 `0.0`，走上面这条「在面内」分支，得到 `Face + hitBoundary=false`，终点其实在边界上却没有边界标记；若算出 `-1e-17` 这类负值，则走 2.2 的跨越分支，得到 `Edge + hitBoundary=true`。
   **实测**：同一平面网格，起始 `(0.25,0.25,0.5)`、位移 `(0.25,0.25,-0.5)`（恰好落到边界边 e0 中点）→ `Face`、`bary=(0.5,0.5,0.0)`、`hitBoundary=false`；而用「重心 → 边界边中点、长度为欧氏范数」的方式瞄准同一条边 → `Edge t=0.5`、`hitBoundary=true`。二者只差浮点舍入的方向。

> 实践含义：**任何「恰好落在边/顶点上」的终点都不能靠 `endPoint.type` 可靠区分。**`hitBoundary==false` 不能推出「终点不在边界上」；要判断终点是否贴边，必须另外检查 `endPoint.type==Face && min(faceCoords)<=0`。

### 2.2 跨越一条边

若终点不在闭三角形内，就在三条边里挑交点（`:214-245`）：

- 只考虑 `edgeIsHittable[i]` 且 `vecBary[i] < 0` 的边（`:233`）；
- 取 `tRayThisRaw = -startPoint[i] / vecBary[i]` 最小者（`:219`、`:239-244`）；
- `tRay` 先被 clamp 到 `[0, 1 - TRACE_EPS_LOOSE]`（`:220`、`:255`），`TRACE_EPS_LOOSE = 1e-9`（`:27`）。上界 `1 - 1e-9` 让交点不会精确落在本面的远端顶点上，也保证剩余长度 `(1-tRay)*len > 0`；但它**不**保证 `tRay > 0`，因此不能由此推出「每步必前进」。
- 交点重心坐标 `endPointOnEdge = startPoint + tRay*vecBary`（`:276`），边参数 `tCross` 由两个分量归一化得到并 clamp 到 `[0,1]`（`:277-283`）。

若一条可跨的边都没选出来（前置条件被破坏），走兜底分支（`:258-273`）：直接以 `SurfacePoint(face, startPoint)` 终止，`traceVectorInHalfedgeLen = 0`。这会给出一个「长度全额、原地未动」的假成功；`errorOnProblem=true` 时改为抛 `std::logic_error`。

### 2.3 边的方向约定与 `tCross`

`SurfacePoint(Edge, tEdge)` 的 `tEdge` 以 `edge.halfedge().vertex()` 为 0（`surface_point.h:25`、`surface_point.ipp:10-12`）。跨边时内部用 `crossHe` 的 `t`，转成规范坐标要经过 `convertTToEdge`（`:87-90`）或在构造 `SurfacePoint(Halfedge,tHalfedge)` 时自动翻转（`surface_point.ipp:11-12`）。

**实测**：4860 个由 `pathPoints` 输出的边点，其 `tEdge` 全部落在 `[0,1]`，且与 `SurfacePoint(edge, tEdge)` 完全一致（不一致 0 个）。所以**对外输出的 `tEdge` 统一是规范 `Edge` 方向**，与内部用哪个 `crossHe` 无关。

**为什么边界边的「来向半边」总是规范半边**：`traceGeodesic` 实际只能跑在 `ManifoldSurfaceMesh` 上——`computeHalfedgeVectorsInVertex` 对 `!usesImplicitTwin()` 直接抛异常（`intrinsic_geometry_interface.cpp:369-371`），而 `traceGeodesic` 无条件 require 它（`:690`）。隐式 twin 网格中 `eHalfedge(iE) = 2*iE`、`heTwinImplicit(iHe) = iHe ^ 1`（`surface_mesh.ipp:58-60`），偶数索引是「真实」半边，边界环半边是奇数索引的隐式 twin 且被标为非 interior（`manifold_surface_mesh.cpp:399-407`）。因此 `boundary edge` 的 `edge.halfedge()` 必为内部半边。**实测**：在 1997 条边界边上（含随机面顺序、随机顶点轮转、单三角形、扇形盘）规范半边为 exterior 的数量为 0。

---

## 3. 穿边时的终点与剩余长度

跨越不会终止追踪（除命中边界/barrier，见 §6）。`traceInFaceBarycentric` 在跨越分支里计算：

```cpp
double lenRemaining = (1.0 - tRay) * vecCartesianLen;                 // :286
Vector2 crossingEdgeVec = vertexCoords[(iOppVertEnd+2)%3]
                        - vertexCoords[(iOppVertEnd+1)%3];           // :287
Vector2 remainingDirInHalfedge = vecCartesianDir / crossingEdgeVec.normalize(); // :288
```

- **剩余长度** `lenRemaining`：本面内按 `tRay` 线性折扣的剩余部分。它就是下一步 `traceInFaceFromEdge` 的 `traceVecInHalfedgeLen`。
- **剩余方向** `remainingDirInHalfedge`：单位复数，坐标系 x 轴是 `crossingEdgeVec`（= `crossHe` 的方向）。**注意这里没有做角度展开/投影**，只是复除法。真正把方向「修正」成合法出射方向的是下一步 `traceInFaceFromEdge`（`:418-446`）：
  - 先对 `traceVecInFaceHalfedgeDir.y = max(y, TRACE_EPS_LOOSE)`（`:420`）强制指向面内；
  - 再转成重心位移并强制 `z >= TRACE_EPS_TIGHT = 1e-12`（`:435`），最后手工把和调成 0（`:437-445`）。

由于这一步是投影而非精确旋转，**每步的实际位移长度与 `lenRemaining` 有微小偏差**；但因为位移向量是用 `vecBaryCanonicalFixedDir * traceVecInHalfedgeLen` 缩放的（`:467-468`），量级仍由 `lenRemaining` 控制。**实测**：请求长度 5.0、穿过 2 条边后命中边界，报告长度与返回路径的三维折线长度逐位相同（差 0）。

---

## 4. 过顶点

### 4.1 从顶点出发

`traceGeodesic_fromVertex`（`:523`）：

1. 绕顶点遍历所有「楔」(wedge)，用 `cross(intervalStart, currVecDir) > 0 && cross(intervalEnd, currVecDir) <= 0` 找出方向所属的楔（`:565-573`）。
2. 若一个都没命中（数值退化），取叉积绝对值最小的边，沿该边方向追踪（`:576-598`），`t` 从 0 出发。
3. 命中楔后：
   ```cpp
   double sum = currVert.isBoundary() ? M_PI : 2. * M_PI;          // :605
   traceDirRelativeToStart = traceDirRelativeToStart.pow(
                                geom.vertexAngleSums[currVert] / sum);  // :606
   ```
   即把输入角度按 `angleSum / (2π)`（内部顶点）或 `angleSum / π`（边界顶点）做**径向缩放**，几何角度 = 输入角度 × 缩放因子。

**实测（平面网格内部顶点 v4，`angleSum=2π` 因子=1）**：输入角与几何角一一对应，只差一个固定旋转（坐标系 x 轴为 `v.halfedge()`）。**实测（锥点，`angleSum=3.486 ≠ 2π`，因子 0.5548）**：输入角被压扁，这正是「把锥面展开成平面」的幂映射。

**实测（边界顶点 v0，`angleSum=π/2`，因子=0.5）**：输入 `[0, π]` 映射到几何 `[0, π/2]`（`π/2` 正好是该顶点两条边界边之一的方向）。输入超出 `[0, π]` 时，第 1 步的楔测试全部失败，落入第 2 步兜底：取叉积绝对值最小的入射半边，从该顶点、沿该半边方向（带 `TRACE_EPS_TIGHT` 的内偏）出发——实测输入 `225°` 得到几何 `0°`、输入 `270°` 得到 `45°`、输入 `360°` 得到 `90°`，都恰好是某条入射边的方向。

所以：**在边界顶点上，只有落在入射楔角度范围内的方向才按「缩放后的角度」解释；超出范围的输入被静默改写为「沿某条入射边出发」，既不报错也不可预测会选哪条。**

### 4.2 途经顶点（起点不是顶点）

这是最容易出问题的情形。追踪器**从不产生 Vertex 类型的中间点**，而是用一串 `t=1`/`t=0` 的 Edge 点绕着顶点「走过去」。机制如下：

- 若精确瞄准顶点且长度恰好耗尽，`isInsideTriangle` 在闭集意义下通过 → 直接以 `Face` 类型终止，`faceCoords=(1,0,0)`（见 §2.1 推论 2）。
- 若越过顶点，下一步的面内终点重心坐标会出现负分量 → 走跨越分支；`tRay` 被 clamp 到 `1 - 1e-9`，产生一个几乎正好在顶点上的 Edge 点；进入相邻面后重复，于是沿顶点扇区逐个楔前进。

**实测（平面，面 0 重心精确瞄准 v4 并多走 0.05）**：
```
path: F0(0.333,0.333,0.333) E1(t=1.000) E8(t=1.000) E14(t=1.000) F7(0.911,0.045,0.045)
三维点: (0.333,0.167) (0.500,0.500) (0.500,0.500) (0.500,0.500) (0.522,0.545)
```
三个 Edge 点在三维上是**同一个点**（顶点 v4），即 `pathPoints` 里存在连续重复的 3D 点。**实测**：该路径 `nPath=5`，连续重复 3D 点 1 处。

**过顶点代价与顶点度相关**。**实测**：在平面扇形网格上从面内精确瞄准中心顶点、多走 0.5，需要跨越的边数：

| 中心顶点度 N | 3 | 4 | 5 | 6 | 8 | 12 | 20 | 40 |
|---|---|---|---|---|---|---|---|---|
| 完整走完所需 `maxIters` | 2 | 2 | 2 | 3 | 4 | 6 | 10 | 20 |

即**穿过一个顶点约需 O(valence/2) 步**。这条直接影响 `maxIters` 的取值（§7）。

**实测**：从内部顶点出发、绕 720 个方向扫描，`pathPoints` 里出现的 Vertex 类型点**只有初始点本身**（720 个），中间点全是 Edge/Face。对 16000 条随机测地线扫描同样如此。

---

## 5. 起点在边上的特殊情形

`traceGeodesic_fromEdge`（`:473`）以「边的切空间」（x 轴 = `edge.halfedge()` 方向，y 轴 = 逆时针 90°）解释输入方向，用 `currVecDir.y >= 0` 决定往哪一侧走（`:486`）：

- `y >= 0`：翻到 `currEdge.halfedge().twin()`，方向取反，`tEdge = 1 - tEdge`（`:487-489`），然后正常进入相邻面。
- `y < 0`：用 `currEdge.halfedge()`；**先检查这条半边是否为边界或 barrier**（`:494-501`），是的话**立即终止**。

**对网格边界边**，因为 `edge.halfedge()` 恒为内部半边（§2.3），`y > 0` 一侧就是相邻面所在侧，`y < 0` 就是曲面外侧。所以这里的行为在几何上是**对的**：朝外走立即停在起点，朝内走正常进入。**实测**（边界边 e0，t=0.1…0.9 共 5 个点）：`y=+0.25` 每次都精确移动 `(0, +0.25)`；`y=-0.25` 每次都停在原处、`endPoint` 仍是同一条边同一个 `t`，`pathPoints` 为 `[E0(t), E0(t)]`。

**问题在长度记账**。立即终止分支构造：

```cpp
result.terminated = true;
result.endPoint = SurfacePoint(currEdge, tEdge);
result.incomingDirToPoint = currVecDir;
// 注意：没有给 result.traceVectorInHalfedgeLen 赋值
```

`TraceSubResult` 是普通结构体，**没有默认成员初始化**（`:103-118`）。主循环随后执行 `result.length -= prevTraceEnd.traceVectorInHalfedgeLen`（`:675`），读到的就是未初始化内存。

**实测**：从边界边中点朝外、请求长度 0.3：
- 普通调用：`length = 0.29999999999999999`（碰巧栈上残留为 0，看起来「正确」，其实是巧合）；
- 调用前用模式值污染栈：`length = -7.6999999999999998e+307`（重复 3 次稳定复现）。

正确语义应当是「在起点即终止，走过长度 0」。**这是一处真实的未初始化读取缺陷**，任何依赖「起点即被边界挡住」时 `length` 值的行为都不可靠。

**对内部 barrier 边**，这里没有「外侧」的概念，`edge.halfedge()` 只是任选的一条内部半边，于是 `y` 的符号只区分「规范半边的哪一侧」，而两侧都是真实的面。结果就是同一 barrier 边、同一 `t`，两侧方向一个停一个走（见 §12 第 2 条）。

另外，这条路径的 `pathPoints` 会是 `[起点, 同一个点]`（`:650-653` 与 `:667-670` 各 push 一次），包含重复点。

---

## 6. 命中边界 / barrier

### 6.1 判据

跨越分支在选择出 `crossHe` 之后检查（`:304`）：

```cpp
if (!crossHe.twin().isInterior() ||
    (traceOptions.barrierEdges && (*traceOptions.barrierEdges)[crossHe.edge()])) {
  result.terminated = true;
  result.endPoint = SurfacePoint(crossHe.edge(), convertTToEdge(crossHe, tCross));
  result.incomingDirToPoint = remainingDirInHalfedge;
  result.traceVectorInHalfedgeLen = lenRemaining;
  return result;
}
```

即「**跨过去的对面不是内部面**」或「**这条边被标记为 barrier**」二者之一。

### 6.2 终点

- 类型恒为 `Edge`；
- 边参数 `tEdge` 已转成规范 `Edge` 方向（`convertTToEdge`）；
- 对真正的网格边界边，`edge.isBoundary()==true`；对 barrier，一般是内部边，`edge.isBoundary()==false`。

**实测（内部 barrier 边）**：`endType=Edge`、`hitBoundary=1`、`endEdge.isBoundary()=false`。这正是 `hitBoundary` 这个名字的误导之处。

### 6.3 `hitBoundary` 的精确含义

`hitBoundary` 只在主循环**正常终止路径**的末尾设置（`:679-681`）：

```cpp
if (prevTraceEnd.endPoint.type == SurfacePointType::Edge) {
  result.hitBoundary = true;
}
```

所以：

> **`hitBoundary == true` ⟺ 追踪因「在边上停下」而终止（边界边或 barrier 边），与边本身是否是网格边界无关。**

反过来，`hitBoundary == false` 有三种可能：(a) 正常耗尽长度（`Face`）；(b) 被 `maxIters` 截断（`Edge`）；(c) 精确落在边界上但被闭集判定吃进了面内（`Face`）。

### 6.4 从 barrier 边上出发

`barrierEdges` 只在「跨越」时检查，但 `traceGeodesic_fromEdge` 对 `y < 0` 的方向会先检查起始边（`:494`）。**实测**：起点在内部 barrier 边上，
- 朝 `y>0` 方向（翻到 twin 侧）：**不触发** barrier，正常走 0.3；
- 朝 `y<0` 方向：**触发** barrier，在起点立即终止。

也就是说「barrier 是否作用于起点」取决于方向落在哪一侧的半边，行为不对称。

---

## 7. 最大步数 `maxIters`

### 7.1 语义

`TraceOptions::maxIters` 默认 `INVALID_IND = SIZE_MAX`（`utilities.h:28`），即**无上限**。

主循环（`:631-665`）：

```cpp
size_t iter = 0;
while (!prevTraceEnd.terminated) {
  if (traceOptions.maxIters != INVALID_IND && iter >= traceOptions.maxIters) {
    result.endPoint   = SurfacePoint(prevTraceEnd.crossHe, prevTraceEnd.tCross);
    result.endingDir  = prevTraceEnd.traceVectorInHalfedgeDir;
    result.length    -= prevTraceEnd.traceVectorInHalfedgeLen;
    if (traceOptions.includePath) result.pathPoints.push_back(result.endPoint);
    return;                      // 提前返回，hitBoundary 保持 false
  }
  ...
  prevTraceEnd = traceInFaceFromEdge(...);
  iter++;
}
```

- `iter` 统计的是**「进入下一个面并走到它的出口」的次数**，即**跨边次数**。
- 起点所在的第一个面走到第一个交点**不在计数内**（那是 `traceGeodesic_from*` 预先完成的）。
- 因此 **`maxIters = k` 表示：第一个交点之后，最多再做 k 步跨越；截断发生在第 `k+1` 个交点处。**

**实测（平面 2×2 网格，请求长度 5，从左边界穿出）**：

| maxIters | 终点 | length |
|---|---|---|
| 0 | 第 1 个交点 e1 | 0.235702260396 |
| 1 | 第 2 个交点 e8 | 0.471404520791 |
| 2 | 第 3 个交点 e13（**自然命中边界**） | 0.942809041582 |
| ≥2 | 同上（已自然终止） | 0.942809041582 |

### 7.2 截断时的终点性质

截断用的是 `prevTraceEnd.crossHe` / `tCross`，而这两个字段**只在「继续」分支被写入**（`:316-321`），即该边必然满足 `crossHe.twin().isInterior() && !barrier`。

> **不变量：`maxIters` 截断的终点一定是「内部、非 barrier」边，绝不会落在边界边或 barrier 边上。**

**实测**（换用 `SurfacePoint` + `Vector2` 重载、方向 `(-1,0)`，因此交点编号与上表不同）：`maxIters=0` 截断于 `edge=2`，`isBoundary=false`、`hitBoundary=false`，请求长度 5.0 只走了 0.1667；`maxIters=1` 起自然命中边界。因此：

> **`endPoint.type == Edge && hitBoundary == false` ⟺ 被 `maxIters` 截断。**

### 7.3 默认无上限时的终止性

默认 `maxIters = SIZE_MAX`，循环只靠几何条件退出。`tRay` 的 clamp 只给上界（`1 - 1e-9`），没有下界，因此**不能**从代码推出「每步都严格消耗长度」，也就**没有**有限步终止的保证。最坏情形是路径反复贴着顶点或落在极扁三角形上。

**实测**：即使瞄准偏移小到 `1e-10`，跨越数也只等于正常过顶点的 `O(valence/2)`，未观察到步数爆炸。但这只是经验观察，不是保证。凡是需要「一定终止」的调用方，都应显式设 `maxIters`。

调用方的实际做法是显式设上限：

- `signpost_intrinsic_triangulation.cpp:77`、`:129`：`options.maxIters = mesh.nFaces() * 10;`
- `intrinsic_triangulation.cpp:269-273`、`integer_coordinates_intrinsic_triangulation.cpp:563-567`：只设 `barrierEdges`，**不设 `maxIters`**。

考虑到过顶点的 `O(valence/2)` 代价，`nFaces * 10` 对常规网格足够；但对高亏格/高 valence 的病态网格并非严格安全。

---

## 8. 剩余长度与 `length` 记账

`result.length` 的生命周期：

1. 入口置为请求长度：`result.length = norm(traceVec)`（`:719`），或重心版 `norm(traceVectorCartesian)`（`:794`）。
2. 单面内终止：`traceVectorInHalfedgeLen = 0`（`:205`、`:271`），故不减。
3. 命中边界/barrier：`traceVectorInHalfedgeLen = lenRemaining`（`:310`），故 `length = 请求长度 - 剩余 = 实际走过多长`。
4. `maxIters` 截断：同 3（`:639`）。
5. 主循环正常终止：`result.length -= prevTraceEnd.traceVectorInHalfedgeLen`（`:675`）。

**实测**：请求 5.0，穿过 2 条边后命中边界，`length = 0.942809041582`，恰等于返回 `pathPoints` 的三维折线长度（差 < 1e-15）。对所有 `maxIters=0/1/2/3` 的组合，报告长度与折线长度逐位一致。

**不变量（实测）**：`length ≤ norm(traceVec)`，最大超出 `1.1e-16`（纯舍入）。截断/边界停止时 `length` 严格更小；**唯一例外是 §5 的未初始化缺陷，此时 `length` 可以是任意垃圾值或负数**。

**特例：零向量。** `traceVec.norm2() == 0` 时提前返回（`:703-717` / `:772-785`）：`length = 0`、`endingDir = (0,0)`，但 **`endPoint` 保持默认构造值**（`type=Vertex`、`vertex` 为空句柄），**不等于起点**；`pathPoints` 里只含起点（当 `includePath=true`）。

---

## 9. 路径完成状态的判据

`TraceGeodesicResult` 的字段（`trace_geodesic.h:12-20`）不直接给出「枚举状态」，需要组合判断。推荐如下判据（每一条分支都有实测覆盖，见 §11）：

```
if (length == 0 && endingDir == 0 && pathPoints.size() <= 1)
    → 零向量提前返回；endPoint 无效，不可使用

else if (endPoint.type == Face)
    if (min(faceCoords) > 0)
        → 正常完成（耗尽长度）；hitBoundary 必为 false；核对 length == 请求长度
    else
        → 精确落在边或顶点上（几何上贴边/在顶点，但按面点返回；
           hitBoundary 可信性见 §2.1，可能为 false）

else if (endPoint.type == Edge)
    if (hitBoundary)
        → 因边而停止：endPoint.edge.isBoundary() ? 网格边界 : barrier 边
           length < 请求长度，等于已走长度
           ⚠ 若起点就在边界/barrier 边上且朝外，length 可能未初始化（§5），不可信
    else
        → 被 maxIters 截断；终点必为内部非 barrier 边
           length < 请求长度，等于已走长度

else if (endPoint.type == Vertex)
    → 不可能是正常追踪的产物（中间点从不产生 Vertex）；
      只会出现在零向量 early-out，或调用方自行构造的结果中
```

补充：

- `hasPath` 只是 `traceOptions.includePath` 的拷贝（`:695`、`:761`），**不表示追踪成功**。
- `pathPoints` 的首元素恒为起点（`includePath` 时）；末元素为 `endPoint`（`:667-670`）。被 `maxIters` 截断时末元素也 push 了（`:641-643`），所以不会少一个点。
- `pathPoints` 中可能存在**连续重复的三维点**（过顶点，§4.2）和**重复的 `(edge,tEdge)`**（从边界边立即终止，§5）。
- `endingDir` 在**正常终止**（面内、边/barrier 停止、`maxIters` 截断）时是单位向量（实测 `|endingDir|=1.000000`），表示**入射方向**，不含长度信息。三个例外：零向量 early-out 时为 `(0,0)`；退化为零长度边时可能是 `inf/NaN`（`:289-300`）；`trimTraceResult` 之后被置 NaN。

### 9.1 `endingDir` 的坐标系（与文档不一致）

头文件称其为「the incoming direction to the final point, in its tangent space」（`trace_geodesic.h:15`），文档进一步规定边的切空间 x 轴为 `edge.halfedge()` 方向。但实现里：

- 命中边界/barrier 时赋值为 `remainingDirInHalfedge`（`:309`），其 x 轴是 **`crossHe`（被跨过的半边）的方向**（`:287-288`）；
- `maxIters` 截断时同样是 `traceVectorInHalfedgeDir`（`:638`），同一坐标系。

`crossHe` 与 `edge.halfedge()` 未必同向。**实测**：对 8 条内部 barrier 边、两侧各测一次，共 16 个样本。跨过的是该边的规范半边（`crossHe == edge.halfedge()`）时，`endingDir` 与规范坐标系一致（误差 ~1e-16）；跨过的是 twin 半边时，`endingDir` 与规范坐标系**相差 180°**（误差 0，与翻转坐标系一致）。即 8 次一致、8 次翻转。

另外，**真正的网格边界边上的停止总是报告为与规范坐标系一致**：对边界边，`edge.halfedge()` 是唯一一条内部半边（§2.3 已证，实测 8/8），所以「来向半边」恒等于规范半边。这也是为什么只看边界情形不会暴露这个问题。

> 因此：**边上的 `endingDir` 用的是「来向半边」坐标系，而不是规范 `Edge` 坐标系**。若需要规范方向，必须自行判断跨过的是哪一侧再决定是否取反。

面内终止时 `endingDir = vecCartesianDir`（`:204`），即面切空间（x 轴 = `face.halfedge()` 方向），与文档一致。**实测**：面内终止 `|endingDir|=1`。

### 9.2 `trimTraceResult`

`trace_geodesic.cpp:811`。用途：当已知追踪「应当」终于某个目标顶点附近时，把尾部的碎点去掉并直接补到该顶点。

行为：

- 从尾部依次弹出：与目标顶点相邻的 Edge 点、所有 Face 点、所有 Vertex 点（`:813-841`），直到遇到一个「不属于上述三类」的点；
- 每弹一个就把 `endingDir = Vector2::undefined()`（NaN）（`:821` 等）；
- 弹出后检查新的尾点是否就是目标顶点、或在目标顶点的 1-ring 内（`:849-875`），返回成功与否；
- **不修改 `endPoint`，也不修改 `length`**；失败时（路径根本没进入目标 1-ring）调用方会重跑一次追踪（`signpost_intrinsic_triangulation.cpp:88-91`、`:139-143`）。

**实测**：对一条 `F E E E F` 的路径 trim 到 v4 → 成功、`pathPoints` 只剩 1 个点、`endingDir` 变为 NaN、`endPoint` 仍是 `Face` 类型（未更新）。trim 到远处顶点 v2 → 失败、`pathPoints` 仍 4 个点。

---

## 10. 曲面内蕴直行 vs 三维直线

这是本题的核心概念区分，实现与实测都支持同一结论。

### 10.1 内蕴直行的定义方式

追踪器在每一步做的是：在**当前面**的二维坐标里沿直线前进，跨边时用 `tCross` 参数定位，然后用「剩余方向 / 边方向」的复除法把方向带进下一个面（`:287-288`），再投影到合法扇形（`:418-446`）。这等价于**把共边三角形依次展开到同一平面后走直线**，即经典的内蕴测地线（只在边长度意义下定义）。

它**没有**做「保持三维方向不变」或「贴合三维表面」的任何计算：代码里没有任何顶点坐标的参与（§1.1）。

### 10.2 与三维长度的关系

因为每个三角形在三维里是**等距嵌入**的（边长由三维坐标算出），所以：

> 报告长度 = 内蕴弧长 = 各段三维直线段长度之和（三维折线长度）。

但**不等于**「三维空间中起点到终点的直线距离（弦长）」。

**实测（平面网格）**：报告长度 0.723802952255，三维折线长度 0.723802952255，差 `-4.4e-16`；因为面共面，折线本身就是三维直线。

**实测（90° 折面条带）**：报告/内蕴长度 `0.900000000`，三维折线长度 `0.900000000`（差 < 1e-15），而**三维弦长只有 0.757674**；每个内部节点的三维转角达 **76°–90°**。也就是说：**内蕴直线在三维里是折线，每个棱处都有明显折角，折角大小由二面角决定**。

**实测（锥点，过顶点）**：让路径精确穿过一个非平坦内部顶点，
- 平面（K=0，`angleSum=2π`）：入向/出向三维方向夹角 **0.0000°**；
- 帐篷（正曲率，`angleSum=3.486`，K=2.797）：夹角 **163.77°**；
- 鞍点（负曲率，`angleSum=5.217`，K=1.066）：夹角 **78.24°**。

即**在锥点处三维方向会突变，突变量由该点的角亏决定**；而在内蕴展开平面里路径是直的。这正是内蕴直行与三维直线的本质差异。

### 10.3 顶点切空间的「展开」

§4.1 的 `pow(angleSum / 2π)` 就是这件事的公式化：顶点处的切空间是一个总角为 `2π`（边界点为 `π`）的抽象锥，输入方向按 `angleSum/2π` 做径向缩放后映射到真实锥面。

**实测（平面内部顶点 v4，`angleSum=2π`）**：因子恰为 1，输入角与几何角严格只差一个固定偏移（`actualAngle = thetaIn + 0.7854`，误差为 0），即恒等映射。

**实测（锥点，`angleSum=3.486`，因子 0.5548）**：以 45° 为间隔输入 8 个方向，相邻两条出射方向三维夹角为 `0.435751 / 0.424752` 交替，输入间隔为 `0.785398`。**同一面内的相邻方向比值严格等于 `0.554815 = angleSum/2π`**（`0.435751 / 0.785398 = 0.554815`）。跨面时三维夹角小于展开角（因为三维里两个不同面的射线之间的欧氏夹角不承载锥面度规），所以不能直接用三维夹角去验证全部 8 步；但**同一面内的比值精确复现缩放因子**，且 8 步的总和 `8 × 平均` 与 `2π × 0.5548` 一致。

注意这个映射的方向：它把**输入圆的 `2π`** 压到 `angleSum`，所以**输入角度在执行时被缩小**，而不是把锥面「摊开成大圆」。对退化网格（`angleSum` 接近 0）这个因子会趋于 0，所有输入方向都塌缩到同一方向；`angleSum > 2π`（鞍点、负曲率）时因子 > 1，输入角被放大——但输入的 `2π` 周期仍然只覆盖 `angleSum` 的实际角度，即**输入方向空间会多圈覆盖锥面**。

### 10.4 何时三维路径恰好也直

只有两种情况：

1. 全程在同一个平面内（单一面，或所有经过的面共面）——此时内蕴=三维；
2. 恰好沿一条直的三维线走，且该线在每一步都落在面内——一般不会自动发生。

**注意**：`SurfacePoint` 的 `interpolate(vertexPositions)` 给出的是三维位置。如果调用方把 `pathPoints` 直接当三维折线用（例如画线、算三维长度、做碰撞），得到的是**正确的三维折线**；但如果用 `length` 去和三维弦长比较，会得到系统性偏大的结果。

**另一个容易混淆的点**：`traceGeodesic` 与 `MeshRayTracer::trace`（`include/geometrycentral/surface/mesh_ray_tracer.h`）是两套完全不同的语义。后者是三维直线–三角形求交（nanort BVH，返回 `t` / `faceIndex` / `faceCoords`），会穿透折痕；前者是内蕴测地线，会在折痕处拐弯。要「三维直线」必须用后者。

---

## 11. 实测证据汇总

测试网格（全部由 `makeManifoldSurfaceMeshAndGeometry` 构造）：

- **平面 2×2 网格**：9 顶点 / 16 边 / 8 面，1 条边界环，中心 v4 为内部顶点，其余 8 个在边界上。
- **扇形网格**：中心点 + N 个环点，所有面以中心为顶点，用于测过顶点的 valence 依赖。
- **80° 帐篷 / 浅鞍**：平面网格把中心顶点抬升 ∓z，得到正/负曲率锥点。
- **90° 折面条带**：两个矩形共一条棱折起。
- **八面体**：闭合曲面（无边界）。

| 编号 | 检验内容 | 观测结果 |
|---|---|---|
| Z1 | 是否依赖三维坐标 | `VertexPositionGeometry` 与 `EdgeLengthGeometry` 结果逐位相同（差异 0） |
| §2.1 | 精确命中内部边 | 16/16 报 `Face`，`hitBoundary=false` |
| §2.1 / §1.2.1 | 精确命中内部顶点 | 报 `Face`，`faceCoords=(1,0,0)`；再以该点作起点时，与 `Vertex` 类型的同一点方向解释完全不同（12/12 不一致） |
| §2.1 | 精确命中边界边 | 两种结果都出现过：`Face+bary=(0.5,0.5,0)`+`hitBoundary=false`，以及 `Edge t=0.5`+`hitBoundary=true` |
| §3 | 长度守恒 | 报告长度 == 三维折线长度（差 < 1e-15） |
| §4.2 | 过顶点 | Edge 点 `t=1`/`t=0` 串，三维上同一点；无 Vertex 类型中间点 |
| §4.2 | 过顶点步数 | 度 N=40 → 需 20 步；N=12 → 6 步；N=6 → 3 步 |
| §5 | 起点在边界边、朝外 | `length` 读到未初始化值：干净栈 0.3，污染栈 `-7.7e307` |
| §5 | 起点在边界边、朝内 | 5 个 t 值下位移都精确为 `(0,+0.25)` |
| §6.2 | 内部 barrier 边停止 | `Edge`、`hitBoundary=1`、`isBoundary=false` |
| §5 / §6.4 / §12 | 起点在内部 barrier 边 | `y=+0.3` 走满 0.3 且 `hitBoundary=false`；`y=-0.3` 立即停在原边、`hitBoundary=true` |
| §7.1 | `maxIters=k` | 停在第 k+1 个交点；`maxIters=0` 停在第一个交点 |
| §7.2 | 截断终点性质 | 必为内部非 barrier 边；`hitBoundary=false` |
| §9 | 零向量 | `length=0`、`endingDir=0`、`endPoint` 为默认无效 `Vertex` |
| §9.1 | 边上的 `endingDir` | 8/16 与规范坐标系一致（规范半边侧），8/16 相差 180°（twin 侧）；边界边停止恒一致 |
| §10.2 | 内蕴 vs 三维 | 折面中折线长=内蕴长，弦长更短；棱处转角 76°–90° |
| §10.2 | 过锥点 | 三维转角：平面 0°，帐篷 163.77°，鞍点 78.24° |

---

## 12. 已知不一致与陷阱清单

1. **未初始化读取（真 bug）**：`traceGeodesic_fromEdge` 的立即终止分支（`:493-501`）不写 `traceVectorInHalfedgeLen`，而 `TraceSubResult` 无默认成员初始化（`:103-118`）。主循环 `:675` 据此计算 `length` → 读到栈上垃圾值。**实测**：从边界边朝外（`y<0`）、请求 0.3，干净栈上 `length=0.3`（巧合），污染栈后稳定为 `-7.7e307`；正确语义应为「起点即终止，走过 0」。

   这也意味着**只有从边界边上朝外出发这一类调用**才可能踩到它；朝内（`y>=0`）走的是另一分支，不经过这里。**实测**：同一起点 `y=+0.25` 正常走 0.25 进入内部，`y=-0.25` 立即停在原边。

2. **`barrierEdges` 对起点不对称**：起点落在 barrier 边上时，`y<0` 分支会先检查该边并立即终止（`:494`），`y>=0` 分支不做这个检查（`:486-489`），于是同一 barrier 边、同一 `t`，两侧方向一个停一个走。**实测**（内部 barrier 边 e1）：`y=+0.3` → `Face`、`hitBoundary=false`、走满 0.3；`y=-0.3` → `Edge`、`hitBoundary=true`、终点仍在 e1 上。若语义是「barrier 边不可穿越」，这里从起点所在边出发本就不构成「穿越」，两侧都不该停；若语义是「barrier 是墙」，则两侧都该停。当前实现两边都不满足，只是恰好各占一半。

3. **`hitBoundary` 名不符实**：它表示「停在边上」，barrier 停止也会置 true（`:679-681`）；且它**不**表示终点在网格边界上。
4. **`hitBoundary=false` 的歧义**：可能是正常完成、`maxIters` 截断、或「精确落在边界但被闭集判据吃进面内」三者之一。
5. **精确落边/落点无类型标记，且类型与几何位置可能不一致**：`isInsideTriangle` 是闭集判据（`>= 0`），导致「几何上在边上/顶点上」的点可能报成 `Face`。调用方若需要区分，必须检查 `min(faceCoords) <= 0`。更麻烦的是，`SurfacePointType` 同时决定切空间（§1.2.1）：把一个顶点写成 `SurfacePoint(face, bary)` 会让输入方向按面坐标解释，结果完全不同且不报错。
6. **过顶点产生重复点**：`pathPoints` 中会出现 3D 坐标完全相同的连续点，以及同一条边的重复 `(edge,tEdge)`。做折线渲染/长度累加时需要去重，否则会得到零长度段（长度不受影响，但拓扑/邻接判断会出错）。
7. **顶点切线空间的角度缩放**：内部顶点按 `angleSum/2π`、边界顶点按 `angleSum/π` 缩放（`:605-606`）。对边界顶点，超出入射楔范围的方向会退化为「沿某条入射边追踪」（`:592-598`），而不是报错。**实测**：边界顶点 v0 上输入 `225°`→几何 `0°`、输入 `270°`→`45°`、输入 `360°`→`90°`。
8. **`endingDir` 的坐标系**：边上终止用「来向半边 `crossHe`」坐标系，与文档/规范 `Edge` 坐标系可能相差 180°（`:287-288`、`:309`、`:638`）。这个字段无法从结果里反查 `crossHe` 是哪条半边，所以**调用方无法自行纠正**，只能接受「方向正确、符号不明」。
9. **`trimTraceResult` 不更新 `endPoint` 与 `length`**，并把 `endingDir` 置为 NaN（`trace_geodesic.h:51-53` 亦有说明）。
10. **零向量 early-out 的 `endPoint` 无效**，不等于起点（`:703-717`、`:772-785`）。
11. **闭合曲面上默认无步数上限**：只有耗尽长度才停；默认 `maxIters = SIZE_MAX`，无有限终止保证。调用方应显式设 `maxIters`（参考 `signpost_intrinsic_triangulation.cpp:77` 的 `nFaces*10`）。

---

## 13. 使用建议（按目标反推）

- **要「走固定长度、拿终点」**：直接用 `endPoint`；用 `length` 与 `norm(traceVec)` 是否相等判断是否走满；`length` 更小时再分 `hitBoundary` 与 `endPoint.type` 两种情况。
- **要判断是否撞墙**：不要只看 `hitBoundary`；对 `Face` 终点还要检查 `min(faceCoords) <= 0`，因为精确贴边可能被吃进面内。
- **要区分边界与 barrier**：`hitBoundary && !endPoint.edge.isBoundary()` 即 barrier。
- **要 `maxIters`**：至少取 `nFaces*10`；若网格有高 valence 顶点（过顶点代价 `O(valence/2)`）或极扁三角形，应更大。
- **要用 `endingDir`**：面内终止时它在面切空间、可信；在边上终止时它可能相对规范 `Edge` 坐标系翻转 180°，而结果里没有信息能判断是否翻转，因此**不要用它来做方向敏感的决策**；`trimTraceResult` 之后该字段为 NaN，不可用。
- **要 3D 路径**：用 `pathPoints[i].interpolate(geometry.vertexPositions)`，并自行合并连续重复点。
- **构造起点**：让 `SurfacePointType` 与几何位置一致（顶点用 `SurfacePoint(v)`，边点用 `SurfacePoint(e,t)`，面内点用 `SurfacePoint(f,bary)`）；不要用「恰好落在顶点上的面点」代替顶点，两者的切空间不同（§1.2.1）。
- **要三维直线**（例如射线/可见性）：不要用 `traceGeodesic`，用 `MeshRayTracer`（`include/geometrycentral/surface/mesh_ray_tracer.h`），它是真正的三维直线–三角形求交，与内蕴测地线是两套语义。
- **需要「精确终点在某顶点」**：`traceGeodesic` 不会返回 `Vertex` 类型；用 `trimTraceResult` 把尾点吸附到目标顶点，并注意它不更新 `endPoint`。
