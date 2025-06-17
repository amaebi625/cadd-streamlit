import os
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import QED, Descriptors, Crippen, Lipinski
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import DataStructs

plt.rcParams['font.sans-serif'] = ['SimHei']  # 显示中文标签
plt.rcParams['axes.unicode_minus'] = False   # 正确显示负号

st.set_page_config(page_title="💊 药物性评估", page_icon="💊")
st.title("💊 药物性评估平台")
st.markdown("本模块用于评估上传分子的基本药物性，包括 **QED 打分**、**Lipinski Rule 检查**、**理化性质计算**等，支持结构图展示与指标可视化。")

# ========== 1. 数据输入 ========== #
st.header("📥 分子数据输入")

input_mode = st.radio("请选择输入方式", ["上传 CSV 文件", "手动输入 SMILES 列表"], horizontal=True)
df = None

if input_mode == "上传 CSV 文件":
    upload = st.file_uploader("上传 CSV 文件（需含 'SMILES' 列）", type=["csv"])
    if upload:
        df = pd.read_csv(upload)
        if 'SMILES' not in df.columns:
            st.error("❌ CSV 文件中必须包含 'SMILES' 列")
            df = None
        else:
            st.success(f"✅ 成功加载 {len(df)} 个分子")
elif input_mode == "手动输入 SMILES 列表":
    smiles_input = st.text_area("请输入 SMILES（每行一个）", height=200,
                                placeholder="CCO\nC1=CC=CC=C1\nO=C(N)Cc1ccccc1")
    if smiles_input.strip():
        smiles_list = smiles_input.strip().splitlines()
        df = pd.DataFrame({"SMILES": smiles_list})
        st.success(f"✅ 成功输入 {len(df)} 个分子")

# ========== 2. 计算指标 ========== #

if df is not None:
    st.header("🧪 指标计算与药物性评估")
    smiles = df["SMILES"].tolist()

    results = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            qed_val = QED.qed(mol)
            mw = Descriptors.MolWt(mol)
            logp = Crippen.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            hba = Lipinski.NumHAcceptors(mol)
            rot = Lipinski.NumRotatableBonds(mol)

            lipinski = int((mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10))

            results.append({
                "SMILES": smi,
                "QED": round(qed_val, 3),
                "MW": round(mw, 1),
                "LogP": round(logp, 2),
                "HBD": hbd,
                "HBA": hba,
                "Rotatable": rot,
                "Lipinski_Pass": "✔️" if lipinski else "❌"
            })
        else:
            results.append({
                "SMILES": smi,
                "QED": None,
                "MW": None,
                "LogP": None,
                "HBD": None,
                "HBA": None,
                "Rotatable": None,
                "Lipinski_Pass": "无效分子 ❌"
            })

    result_df = pd.DataFrame(results)
    st.dataframe(result_df)

    # ========== 2.5 QED筛选功能 ========== #
    st.subheader("🔍 QED筛选")
  
    if len(result_df) == 0:
        st.warning("⚠️ 当前没有可用于筛选的分子，请先完成药物性评价步骤。")
        st.stop()
    
    filter_mode = st.radio("选择筛选方式", ["按QED阈值", "保留Top-N高分子"])
    
    if filter_mode == "按QED阈值":
        qed_threshold = st.slider("QED分数阈值", 0.0, 1.0, 0.6)
        filtered_df = result_df[result_df["QED"] >= qed_threshold]
    else:
        max_n = len(result_df)
        default_n = min(10, max_n)
        top_n = st.number_input(
            "保留前 N 个QED最高分子",
            min_value=1,
            max_value=max_n,
            value=default_n,
            step=1
        )
        filtered_df = result_df.sort_values("QED", ascending=False).head(top_n)
    
    st.success(f"✅ 筛选后剩余 {len(filtered_df)} 个分子")
    st.dataframe(filtered_df)
    
    # ✅ 保存筛选结果到 session_state，用于后续活性预测
    if st.button("✔️ 使用当前筛选结果进行性质预测"):
        st.session_state["selected_for_activity"] = filtered_df.copy()
        st.success(f"✅ 已保存 {len(filtered_df)} 个分子用于性质预测，可前往下一页查看")

    # ========== 3. 可视化雷达图 ========== #
    st.subheader("📊 雷达图展示(最多5个分子)")
    selected = filtered_df.dropna().head(5)

    def normalize_properties(row):
        return {
            'MW': row['MW'] / 500,
            'LogP': row['LogP'] / 5,
            'HBD': row['HBD'] / 5,
            'HBA': row['HBA'] / 10,
            'Rotatable': row['Rotatable'] / 10,
            'QED': row['QED']
        }

    def plot_multi_radar(df, title="Multi-Sample Radar Chart"):
        normed_rows = [normalize_properties(row) for _, row in df.iterrows()]
        radar_df = pd.DataFrame(normed_rows)
        labels = radar_df.columns.tolist()
        angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False).tolist()
        angles += angles[:1]

        fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

        for i, row in radar_df.iterrows():
            values = row.tolist()
            values += values[:1]
            ax.plot(angles, values, label=f"Sample {i+1}")
            ax.fill(angles, values, alpha=0.1)

        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.set_thetagrids(np.degrees(angles[:-1]), labels, fontsize=10)
        ax.set_title(title)
        ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
        return fig

    if selected.empty:
        st.info("⚠️ 筛选结果为空，暂无可视化分子")
    else:
        fig = plot_multi_radar(selected)
        st.pyplot(fig)

    # ========== 4. 下载按钮 ========== #
    st.download_button(
        label="📥 下载评估结果 CSV",
        data=filtered_df.to_csv(index=False).encode("utf-8"),
        file_name="drug_eval_result.csv",
        mime="text/csv"
    )
