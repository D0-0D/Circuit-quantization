import streamlit as st
import sympy as sp
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# 导入你封装好的引擎模块
from build_circuit_graph import Component, Circuit, MutualInductance

# ==================== 页面配置 ====================
st.set_page_config(page_title="超导电路量子化工具", layout="wide")
st.title("⚛️ 超导电路量子化与哈密顿量生成器")

# 初始化 session_state，用于存放默认的电路表格数据
if 'components_df' not in st.session_state:
    st.session_state.components_df = pd.DataFrame({
        'u': [0, 1, 0, 0, 1],
        'v': [2, 2, 1, 2, 2],
        'type': ['JJ', 'JJ', 'L', 'C', 'C'],
        'value': ['EJ1', 'EJ2', 'L', 'C1', 'C2']
    })

# ==================== 第 1 步：创建电路对象 ====================
st.header("1. 定义电路元件")
st.markdown("请在下方表格中添加、删除或修改元件参数（点击表格即可编辑）：")

# 使用交互式表格编辑元件
edited_df = st.data_editor(
    st.session_state.components_df, 
    num_rows="dynamic",
    column_config={
        "u": st.column_config.NumberColumn("节点 u", step=1),
        "v": st.column_config.NumberColumn("节点 v", step=1),
        "type": st.column_config.SelectboxColumn("元件类型", options=["C", "L", "JJ"], required=True),
        "value": st.column_config.TextColumn("参数值 (符号/数值)", required=True)
    },
    use_container_width=True
)

# 尝试根据表格数据实例化电路
circuit = None
try:
    my_components = []
    for index, row in edited_df.iterrows():
        # 忽略空行
        if pd.isna(row['u']) or pd.isna(row['v']) or pd.isna(row['type']):
            continue
        my_components.append(Component(int(row['u']), int(row['v']), row['type'], str(row['value'])))
    
    if my_components:
        circuit = Circuit(my_components)
        st.success(f"✅ 电路已成功构建！包含 {circuit.nt} 个树枝，{circuit.num_loops} 个基本回路。")
    else:
        st.warning("请添加至少一个元件。")
except Exception as e:
    st.error(f"构建电路时发生错误: {e}")


# 如果电路实例化成功，进行后续步骤
if circuit:
    col1, col2 = st.columns([1, 1])
    
    with col1:
        # ==================== 第 2 步：查看基本回路 ====================
        st.header("2. 查看基本回路")
        
        # 绘制电路图 (使用 NetworkX)
        fig, ax = plt.subplots(figsize=(4, 3))
        # 提取用于绘图的图结构
        G_draw = nx.MultiGraph()
        for k, info in circuit.edge_map.items():
            G_draw.add_edge(info['u'], info['v'], label=f"{info['type']}({info['value']})")
        
        pos = nx.spring_layout(G_draw)
        nx.draw(G_draw, pos, ax=ax, with_labels=True, node_color='lightblue', node_size=500, font_weight='bold')
        
        # 处理多重图的边标签绘制
        edge_labels = {}
        for u, v, data in G_draw.edges(data=True):
            # 简单处理：如果是多重边，标签拼在一起
            if (u, v) in edge_labels:
                edge_labels[(u, v)] += f", {data['label']}"
            else:
                edge_labels[(u, v)] = data['label']
                
        nx.draw_networkx_edge_labels(G_draw, pos, edge_labels=edge_labels, ax=ax, font_size=8)
        st.pyplot(fig)

        # 打印回路信息
        with st.expander("查看所有基本回路详情", expanded=True):
            for i, loop in enumerate(circuit.loops):
                chord = loop['chord_info']
                path_str = ' -> '.join(map(str, loop['path_nodes'])) + f" -> {loop['path_nodes'][0]}"
                st.markdown(f"**回路 {i}** (由连支 `{chord['type']}: {chord['value']}` 闭合)")
                st.markdown(f"- 路径: `{path_str}`")

    with col2:
        # ==================== 第 3 步：编辑外磁通 ====================
        st.header("3. 编辑外磁通")
        st.markdown("为每个基本回路分配外磁通（输入 `0` 表示无外磁通）：")
        
        # 动态生成外磁通输入框，并在 circuit 中进行设置
        if circuit.num_loops == 0:
            st.info("当前电路无闭合回路，不需要设置外磁通。")
        else:
            for i in range(circuit.num_loops):
                # 默认值可以设为 0 或者 Phi_ext_i
                default_val = "0"
                if i == 0: default_val = "Phi_e" # 给第一个回路一个默认的外磁通符号作为演示
                
                flux_input = st.text_input(f"回路 {i} 的外磁通 $\\Phi_{{ext,{i}}}$", value=default_val, key=f"flux_{i}")
                
                # 处理输入：如果输入是纯数字字符串转为 int/float，否则作为符号字符串传入
                try:
                    if '.' in flux_input:
                        flux_val = float(flux_input)
                    else:
                        flux_val = int(flux_input)
                except ValueError:
                    flux_val = flux_input # 保持为字符串，引擎内部会转为 sympy 符号
                    
                circuit.set_external_flux(i, flux_val)

    st.divider()

    # ==================== 第 4 步：计算哈密顿量 ====================
    st.header("4. 计算哈密顿量")
    
    if st.button("🚀 开始计算 (Generate Hamiltonian)", type="primary"):
        with st.spinner("正在进行符号推导与矩阵计算..."):
            try:
                H, info = circuit.hamiltonian()
                
                # 使用两列展示结果，更美观
                res_col1, res_col2 = st.columns(2)
                
                with res_col1:
                    st.subheader("广义坐标 (树枝磁通) $\\Phi_t$")
                    st.latex(sp.latex(info['phi_t']))
                    
                    st.subheader("基本割集矩阵 $F_C$")
                    st.latex(sp.latex(info['F_C']))

                with res_col2:
                    st.subheader("外磁通变量")
                    if info['ext_fluxes']:
                        st.latex(sp.latex(info['ext_fluxes']))
                    else:
                        st.write("无外磁通")
                        
                    st.subheader("全支路磁通 $\\vec{\\Phi}$")
                    st.latex(sp.latex(info['Phi_vec']))
                
                st.subheader("总哈密顿量 $\\hat{H}$")
                st.info("提示：动能项使用了 $M^{-1}$ 表示质量矩阵的逆。")
                
                # 展开后的哈密顿量可能会很长，单独占用一整行
                st.latex(sp.latex(sp.expand(H)))
                
            except Exception as e:
                st.error(f"计算过程中发生数学或图论错误: {e}")