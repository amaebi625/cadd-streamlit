# pages/2_水合自由能预测.py

import streamlit as st
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import matplotlib.pyplot as plt

st.set_page_config(page_title="💧 水合自由能预测", page_icon="💧", layout="wide")
st.title("💧 药物分子的水合自由能预测")

# === 检查 Session 中是否有保留的分子 === #
if "selected_for_activity" not in st.session_state:
    st.warning("⚠️ 请先在【药物性评估】页完成分子筛选。")
    st.stop()

df = st.session_state["selected_for_activity"].copy()
st.success(f"✅ 已载入 {len(df)} 个分子用于预测")

# === 加载模型 === #
try:
    model = joblib.load("config/hydration_rf_model.pkl")
except FileNotFoundError:
    st.error("❌ 未找到训练好的模型文件，请先运行 train_hydration_model.py 脚本")
    st.stop()

# === 分子转 Fingerprint === #
def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None

df["fp"] = df["SMILES"].apply(smiles_to_fp)
df = df[df["fp"].notnull()].copy()

X = np.array(df["fp"].tolist())
pred = model.predict(X)
df["HydrationEnergy (kcal/mol)"] = np.round(pred, 3)

# === 展示预测结果 === #
st.subheader("📈 预测结果展示")
st.dataframe(df[["SMILES", "HydrationEnergy (kcal/mol)"]])

# === 可视化 === #
st.subheader("🎯 水合自由能分布图")
fig, ax = plt.subplots()
ax.hist(df["HydrationEnergy (kcal/mol)"], bins=15, color="#1f77b4", edgecolor="black")
ax.set_xlabel("Predicted Hydration Free Energy (kcal/mol)")
ax.set_ylabel("Count")
st.pyplot(fig)

# === 下载结果 === #
st.download_button(
    label="📥 下载预测结果",
    data=df[["SMILES", "HydrationEnergy (kcal/mol)"]].to_csv(index=False).encode("utf-8"),
    file_name="hydration_energy_pred.csv",
    mime="text/csv"
)
