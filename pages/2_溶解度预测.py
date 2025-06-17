# pages/2_水溶解度预测.py
import streamlit as st
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

st.set_page_config(page_title="💧 水溶解度预测", page_icon="💧", layout="wide")
st.title("💧 药物分子的水溶解度预测（ESOL）")

# ===== 读取第一页保留的分子 =====
if "selected_for_activity" not in st.session_state:
    st.warning("⚠️ 请先在【药物性评估】页完成分子筛选。")
    st.stop()

df = st.session_state["selected_for_activity"].copy()
st.success(f"✅ 已载入 {len(df)} 个分子用于溶解度预测")

# ===== 加载模型 =====
try:
    model = joblib.load("config/esol_rf_model.pkl")
except FileNotFoundError:
    st.error("❌ 未找到训练好的模型文件，请先运行 train_esol_model.py 脚本进行训练。")
    st.stop()

# ===== 特征提取函数 =====
def smiles_to_features(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return [
        Descriptors.MolWt(mol),
        Crippen.MolLogP(mol),
        Lipinski.NumHDonors(mol),
        Lipinski.NumHAcceptors(mol),
        Lipinski.NumRotatableBonds(mol),
        Descriptors.TPSA(mol),
    ]

df["features"] = df["SMILES"].apply(smiles_to_features)
df = df.dropna(subset=["features"])
X = np.array(df["features"].tolist())
df["ESOL_predicted"] = model.predict(X)

# ===== 展示预测结果 =====
st.subheader("📋 溶解度预测结果")
st.dataframe(df[["SMILES", "QED", "MW", "LogP", "HBD", "HBA", "Rotatable", "ESOL_predicted"]])

# ===== 图形展示：MW vs 溶解度 =====
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(figsize=(6, 4))
sns.scatterplot(x="MW", y="ESOL_predicted", data=df, hue="QED", palette="viridis", ax=ax)
ax.set_title("分子量与预测溶解度关系图")
ax.set_xlabel("分子量 (MW)")
ax.set_ylabel("预测 Log(溶解度)")
st.pyplot(fig)

# ===== 下载结果 =====
st.download_button(
    label="📥 下载溶解度预测结果 CSV",
    data=df.to_csv(index=False).encode("utf-8"),
    file_name="solubility_prediction.csv",
    mime="text/csv"
)
