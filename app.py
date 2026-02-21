import streamlit as st
import sympy as sp
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from build_circuit_graph import Component, Circuit, MutualInductance

# ==================== 页面配置 ====================
st.set_page_config(page_title="超导电路量子化工具", layout="wide")
st.title("超导电路量子化与哈密顿量生成器")

# ==================== 初始化 session_state ====================
if 'components_df' not in st.session_state:
    st.session_state.components_df = pd.DataFrame({
        'id': [0, 1, 2, 3, 4],
        'u': [0, 1, 0, 0, 1],
        'v': [2, 2, 1, 2, 2],
        'type': ['JJ', 'JJ', 'L', 'C', 'C'],
        'value': ['EJ1', 'EJ2', 'L', 'C1', 'C2']
    })

if 'mutuals_df' not in st.session_state:
    st.session_state.mutuals_df = pd.DataFrame({
        'L1_id': pd.Series(dtype='int'),
        'L2_id': pd.Series(dtype='int'),
        'value': pd.Series(dtype='str')
    })

# ==================== 第 1 步：定义电路元件 ====================
st.header("1. 定义电路元件")
st.info(
    "**接地节点约定**：编号**最大**的节点默认为接地节点。"
    "例如节点编号为 0, 1, 2 时，节点 2 为接地节点。"
)

st.subheader("元件列表")
st.markdown("请在下方表格中添加、删除或修改元件参数（`id` 列用于互感引用，请保持唯一）：")

edited_df = st.data_editor(
    st.session_state.components_df,
    num_rows="dynamic",
    column_config={
        "id": st.column_config.NumberColumn("元件 ID", step=1, help="元件编号，互感中引用此 ID"),
        "u": st.column_config.NumberColumn("节点 u", step=1),
        "v": st.column_config.NumberColumn("节点 v", step=1),
        "type": st.column_config.SelectboxColumn("元件类型", options=["C", "L", "JJ"], required=True),
        "value": st.column_config.TextColumn("参数值 (符号/数值)", required=True),
    },
    width='stretch',
)

# ==================== 互感定义 ====================
st.subheader("互感列表（可选）")
st.markdown("若电路中存在互感耦合，请在下方填写。`L1_id` 和 `L2_id` 对应上方元件表中电感的 **元件 ID**。")

edited_mutuals_df = st.data_editor(
    st.session_state.mutuals_df,
    num_rows="dynamic",
    column_config={
        "L1_id": st.column_config.NumberColumn("电感1 ID", step=1),
        "L2_id": st.column_config.NumberColumn("电感2 ID", step=1),
        "value": st.column_config.TextColumn("互感值 (符号/数值)", required=True),
    },
    width='stretch',
)

# ==================== 构建电路 ====================
circuit = None
try:
    # 解析元件
    my_components = []
    id_to_index = {}  # 用户指定的 id -> 元件在列表中的位置索引
    for idx, (_, row) in enumerate(edited_df.iterrows()):
        if pd.isna(row.get('u')) or pd.isna(row.get('v')) or pd.isna(row.get('type')):
            continue
        comp_id = int(row['id']) if not pd.isna(row.get('id')) else idx
        id_to_index[comp_id] = idx
        my_components.append(
            Component(int(row['u']), int(row['v']), row['type'], str(row['value']))
        )

    # 解析互感
    my_mutuals = []
    for _, row in edited_mutuals_df.iterrows():
        if pd.isna(row.get('L1_id')) or pd.isna(row.get('L2_id')) or pd.isna(row.get('value')):
            continue
        l1_id = int(row['L1_id'])
        l2_id = int(row['L2_id'])
        # 将用户的元件 ID 映射到列表位置索引（引擎按列表顺序编号）
        if l1_id not in id_to_index or l2_id not in id_to_index:
            st.warning(f"互感引用的元件 ID ({l1_id}, {l2_id}) 在元件表中未找到，已跳过。")
            continue
        my_mutuals.append(
            MutualInductance(id_to_index[l1_id], id_to_index[l2_id], str(row['value']))
        )

    if my_components:
        circuit = Circuit(my_components, my_mutuals if my_mutuals else None)
        st.success(f"✅ 电路已成功构建！包含 {circuit.nt} 个树枝，{circuit.num_loops} 个基本回路。")
    else:
        st.warning("请添加至少一个元件。")
except Exception as e:
    st.error(f"构建电路时发生错误: {e}")

# ==================== 后续步骤（电路构建成功后） ====================
if circuit:

    # ==================== 第 2 步：重排后的元件表 ====================
    st.header("2. 重排后的元件表（树枝优先）")
    st.markdown(
        "引擎会按照 **电容 → 约瑟夫逊结 → 电感** 的优先级选择生成树，"
        "并对所有支路重新编号。下表展示重排后的编号，后续矩阵和公式均使用此编号。"
    )
    reordered_rows = []
    for k in sorted(circuit.edge_map.keys()):
        info = circuit.edge_map[k]
        role = "树枝" if k < circuit.nt else "连支"
        reordered_rows.append({
            'Key': k,
            '角色': role,
            '类型': info['type'],
            '节点 u': info['u'],
            '节点 v': info['v'],
            '参数值': str(info['value']),
        })
    reordered_df = pd.DataFrame(reordered_rows)
    st.dataframe(reordered_df, hide_index=True, width='stretch')

    col1, col2 = st.columns([1, 1])

    with col1:
        # ==================== 第 3 步：查看基本回路 ====================
        st.header("3. 查看基本回路")

        # 绘制电路图 (使用 NetworkX)
        fig, ax = plt.subplots(figsize=(4, 3))
        G_draw = nx.MultiGraph()
        for k, info in circuit.edge_map.items():
            G_draw.add_edge(info['u'], info['v'], label=f"{info['type']}({info['value']})")

        pos = nx.spring_layout(G_draw, seed=42)
        nx.draw(
            G_draw, pos, ax=ax, with_labels=True,
            node_color='lightblue', node_size=500, font_weight='bold',
        )

        # 处理多重图的边标签
        edge_labels = {}
        for u, v, data in G_draw.edges(data=True):
            if (u, v) in edge_labels:
                edge_labels[(u, v)] += f", {data['label']}"
            else:
                edge_labels[(u, v)] = data['label']

        nx.draw_networkx_edge_labels(G_draw, pos, edge_labels=edge_labels, ax=ax, font_size=8)
        st.pyplot(fig)

        # 回路详情
        with st.expander("查看所有基本回路详情", expanded=True):
            if circuit.num_loops == 0:
                st.info("当前电路无闭合回路。")
            for i, loop in enumerate(circuit.loops):
                chord = loop['chord_info']
                path_str = ' -> '.join(map(str, loop['path_nodes'])) + f" -> {loop['path_nodes'][0]}"
                st.markdown(f"**回路 {i}** (由连支 Key {loop['chord_key']}: `{chord['type']}: {chord['value']}` 闭合)")
                st.markdown(f"- 路径: `{path_str}`")

    with col2:
        # ==================== 第 4 步：编辑外磁通 ====================
        st.header("4. 编辑外磁通")
        st.markdown("为每个基本回路分配外磁通（输入 `0` 表示无外磁通）：")

        if circuit.num_loops == 0:
            st.info("当前电路无闭合回路，不需要设置外磁通。")
        else:
            for i in range(circuit.num_loops):
                default_val = "0"
                if i == 0:
                    default_val = "Phi_e"

                flux_input = st.text_input(
                    f"回路 {i} 的外磁通 $\\Phi_{{ext,{i}}}$",
                    value=default_val,
                    key=f"flux_{i}",
                )

                # 解析输入
                try:
                    if '.' in flux_input:
                        flux_val = float(flux_input)
                    else:
                        flux_val = int(flux_input)
                except ValueError:
                    flux_val = flux_input

                circuit.set_external_flux(i, flux_val)

    st.divider()

    # ==================== 第 5 步：计算哈密顿量 ====================
    st.header("5. 计算哈密顿量")

    if st.button("开始计算", type="primary"):
        with st.spinner("正在进行符号推导与矩阵计算..."):
            try:
                H, info = circuit.hamiltonian()

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
                        st.latex("0")

                    st.subheader("全支路磁通 $\\vec{\\Phi}$")
                    st.latex(sp.latex(info['Phi_vec']))

                st.subheader("总哈密顿量 $\\hat{H}$")
                # 仅在质量矩阵奇异（无法求逆）时显示提示
                if sp.Symbol("M^{-1}") in H.free_symbols:
                    st.info("提示：质量矩阵奇异，动能项使用了符号 $M^{-1}$ 表示质量矩阵的逆。")
                st.latex(sp.latex(sp.expand(H)))

            except Exception as e:
                st.error(f"计算过程中发生数学或图论错误: {e}")
