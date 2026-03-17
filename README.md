# 超导电路量子化工具

基于图论（生成树、基本割集）对超导量子比特电路进行**符号化量子化**，自动推导哈密顿量。

## 功能

- 支持电容 (C)、电感 (L)、约瑟夫逊结 (JJ) 三类元件
- 支持电感间的**互感耦合**
- 自动选择最小权重生成树（JJ/L 优先进入树枝，电容倾向连支）
- 计算基本割集矩阵 $Q_f$
- **两种外磁通设置方式**（二选一，不可混用）：
  - **方式 A**：通过**物理回路**绑定外磁通，支持多磁通自动叠加
  - **方式 B**：直接指定基本回路磁通，可指定任意个，未指定的默认为 0
- 符号推导完整哈密顿量：动能（电场能）+ 电感势能 + 约瑟夫逊势能
- 支持 Jupyter Notebook 交互式使用

## 快速开始

### 安装依赖

```bash
pip install -r requirements.txt
```

### 作为 Python 模块使用

```python
from build_circuit_graph_rebuild import Circuit

# 创建空电路
circuit = Circuit()

# 添加元件，返回元件 ID
cj1 = circuit.add_component(0, 2, 'JJ', 'EJ1')
cj2 = circuit.add_component(1, 2, 'JJ', 'EJ2')
cl  = circuit.add_component(0, 1, 'L',  'L')
cc1 = circuit.add_component(0, 2, 'C',  'C1')
cc2 = circuit.add_component(1, 2, 'C',  'C2')

# 添加互感（只能施加在电感之间）
# mid = circuit.add_mutual(cl1, cl2, 'M12')

# ---- 外磁通方式 A：物理回路 ----
# 需提供恰好 num_loops 个线性无关的物理回路（含磁通为 0 的）
circuit.add_physical_flux([cc1, cc2, cl], 'Phi_ext')
circuit.add_physical_flux([cj1, cc1], 0, direction=1)
circuit.add_physical_flux([cj2, cc2], 0, direction=1)

# ---- 外磁通方式 B：直接指定基本回路磁通（与方式 A 二选一）----
# circuit.set_external_flux(0, 'Phi_ext')   # 只需设置非零的回路
# circuit.get_external_flux(0)              # 查询
# circuit.clear_external_fluxes()           # 清除后可切换回方式 A

# 查看信息
circuit.print_edges()
circuit.print_loops()
circuit.print_physical_fluxes()

# 计算哈密顿量
H, info = circuit.hamiltonian()
```

## 使用说明

### 元件类型

| 类型 | 说明 | 生成树权重 |
|------|------|-----------|
| `JJ` | 约瑟夫逊结 | 0（最优先进入树枝） |
| `L`  | 电感 | 1 |
| `C`  | 电容 | 2（倾向连支） |

### 约定

- **接地节点**：编号最大的节点为接地节点
- **边方向**：始终从小编号节点指向大编号节点
- **参数值**：支持符号字符串（如 `'EJ1'`）或数值
- **元件 ID**：`add_component` 返回的整数 ID，用于后续引用（互感、磁通、删除）

### 互感

通过元件 ID 添加互感，只能施加在电感 (`L`) 之间：

```python
cl1 = circuit.add_component(0, 1, 'L', 'L1')
cl2 = circuit.add_component(1, 2, 'L', 'L2')
circuit.add_mutual(cl1, cl2, 'M12')
```

### 外磁通

提供两种互斥的设置方式：

#### 方式 A：物理回路磁通

通过元件 ID 列表指定物理回路，顺序决定磁通正方向（右手定则）。需提供恰好 `num_loops` 个线性无关的物理回路（含磁通为 0 的），相同物理回路的磁通自动合并：

```python
circuit.add_physical_flux([cj1, cl, cj2], 'Phi_e')   # 有磁通的回路
circuit.add_physical_flux([cj1, cc1], 0, direction=1)  # 磁通为 0 也要添加
```

#### 方式 B：直接指定基本回路磁通

直接为基本回路设置外磁通值，可指定任意个，未指定的默认为 0：

```python
circuit.set_external_flux(0, 'Phi_ext')  # 回路 0 的外磁通
circuit.set_external_flux(1, 0)          # 可选：显式设为 0
circuit.get_external_flux(0)             # 查询回路 0 的磁通
circuit.fundamental_fluxes               # 查看所有已设置的 {回路编号: 磁通}
circuit.clear_external_fluxes()          # 清除所有（清除后可切换回方式 A）
```

> **注意**：两种方式不可混用。若已使用 `set_external_flux`，需先调用 `clear_external_fluxes()` 才能使用 `add_physical_flux`，反之亦然（需先删除所有物理磁通）。修改元件/图结构时，已设置的基本回路磁通会自动清除。

### 动态修改

元件、互感、磁通均可随时增删，图结构按需自动重建：

```python
circuit.remove_component(cc1)   # 删除元件（关联互感一并删除）
circuit.remove_mutual(mid)      # 删除互感
circuit.remove_physical_flux(0) # 删除磁通
```

## 项目结构

```
Circuit-quantization/
├── build_circuit_graph_rebuild.py  # 核心引擎
├── test_rebuild.ipynb              # 使用示例与测试
├── requirements.txt                # Python 依赖
└── README.md
```

## 算法流程

1. **构建电路图** — 将元件构建为 NetworkX 多重图，计算最小权重生成树，对所有支路重新编号（树枝 `0..nt-1`，连支 `nt..m-1`）
2. **基本割集矩阵** — 对每个树枝，移除后找连通分量，确定哪些连支穿越割集，构建 $Q_f$ 矩阵
3. **参数矩阵** — 构建电容矩阵 $D_C$、电感逆矩阵 $L^+$（含互感）、约瑟夫逊矩阵 $D_J$
4. **外磁通** — 方式 A：将物理回路投影到基本回路空间，解线性方程组；方式 B：直接使用用户指定的值
5. **哈密顿量** — 动能 $\frac{1}{2} q^T M^{-1} q$ + 电感势能 $\frac{1}{2} \Phi^T L^+ \Phi$ + 约瑟夫逊势能 $-\sum E_J \cos(\Phi_k)$
