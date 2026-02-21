# 超导电路量子化工具

基于图论（生成树、基本割集）对超导量子比特电路进行**符号化量子化**，自动推导哈密顿量。

## 功能

- 支持电容 (C)、电感 (L)、约瑟夫逊结 (JJ) 三类元件
- 支持电感间的**互感耦合**
- 自动选择最小权重生成树（电容优先进入树枝）
- 计算基本割集矩阵 $Q_f$
- 为每个基本回路分配外磁通
- 符号推导完整哈密顿量：动能（电场能）+ 电感势能 + 约瑟夫逊势能
- 提供 Streamlit Web 界面，支持交互式编辑和可视化

## 快速开始

### 安装依赖

```bash
pip install -r requirements.txt
```

### 启动 Web 界面

```bash
streamlit run app.py
```

### 作为 Python 模块使用

```python
from build_circuit_graph import Component, Circuit, MutualInductance

# 定义元件
components = [
    Component(0, 2, 'JJ', 'EJ1'),
    Component(1, 2, 'JJ', 'EJ2'),
    Component(0, 1, 'L',  'L'),
    Component(0, 2, 'C',  'C1'),
    Component(1, 2, 'C',  'C2'),
]

# 创建电路
circuit = Circuit(components)

# 查看信息
circuit.print_edges()
circuit.print_loops()

# 设置外磁通
circuit.set_external_flux(0, 'Phi_e')
circuit.set_external_flux(1, 0)

# 计算哈密顿量
H, info = circuit.hamiltonian()
```

## 使用说明

### 元件类型

| 类型 | 说明 | 生成树权重 |
|------|------|-----------|
| `C`  | 电容 | 1（优先进入树枝） |
| `JJ` | 约瑟夫逊结 | 2 |
| `L`  | 电感 | 3 |

### 约定

- **接地节点**：编号最大的节点为接地节点
- **边方向**：始终从小编号节点指向大编号节点
- **参数值**：支持符号字符串（如 `'EJ1'`）或数值

### 互感

定义互感时，`L1_id` 和 `L2_id` 对应元件在输入列表中的索引（从 0 开始）：

```python
mutuals = [
    MutualInductance(L1_id=2, L2_id=5, value='M12'),
]
circuit = Circuit(components, mutuals)
```

## 项目结构

```
Circuit-quantization/
├── build_circuit_graph.py   # 核心引擎：图构建、生成树、割集矩阵、哈密顿量推导
├── app.py                   # Streamlit Web 界面
├── requirements.txt         # Python 依赖
└── README.md
```

## 算法流程

1. **构建电路图** — 将元件列表构建为 NetworkX 多重图，计算最小权重生成树，对所有支路重新编号（树枝 `0..nt-1`，连支 `nt..m-1`）
2. **基本割集矩阵** — 对每个树枝，移除后找连通分量，确定哪些连支穿越割集，构建 $Q_f$ 矩阵
3. **参数矩阵** — 构建电容矩阵 $D_C$、电感逆矩阵 $L^+$（含互感）、约瑟夫逊矩阵 $D_J$
4. **哈密顿量** — 动能 $\frac{1}{2} q^T M^{-1} q$ + 电感势能 $\frac{1}{2} \Phi^T L^+ \Phi$ + 约瑟夫逊势能 $-\sum E_J \cos(\Phi_k)$
