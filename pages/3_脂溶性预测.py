import streamlit as st
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

st.set_page_config(page_title="💧 亲脂性预测 (XGBoost)", page_icon="💧", layout="wide")
st.title("💧 药物分子的亲脂性预测（LogP，XGBoost 模型）")

# ===== 读取上一页筛选的分子 =====
if "selected_for_activity" not in st.session_state:
    st.warning("⚠️ 请先在【药物性评估】页完成分子筛选。")
    st.stop()

df = st.session_state["selected_for_activity"].copy()
st.success(f"✅ 已载入 {len(df)} 个分子用于 LogP 预测")

# ===== 加载 XGBoost 模型 =====
try:
    model = joblib.load("config/logp_xgb_model.pkl")
except FileNotFoundError:
    st.error("❌ 未找到训练好的 XGBoost 模型文件，请先运行 train_logp_model.py 脚本生成 config/logp_xgb_model.pkl。")
    st.stop()

# ===== 分子转指纹函数 =====
def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    arr = np.zeros((2048,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

fps = []
valid_indices = []
for i, smi in enumerate(df["SMILES"]):
    fp = smiles_to_fp(smi)
    if fp is not None:
        fps.append(fp)
        valid_indices.append(i)

if not fps:
    st.error("❌ 没有有效的 SMILES 可用于预测。")
    st.stop()

X = np.array(fps)
logp_preds = model.predict(X)

# ===== 展示预测结果 =====
df_valid = df.iloc[valid_indices].copy()
df_valid["预测 LogP"] = np.round(logp_preds, 2)
st.dataframe(df_valid)

# ===== 可视化 =====
st.subheader("📈 LogP 分布图")
st.bar_chart(df_valid[["预测 LogP"]].reset_index(drop=True))

# ===== 下载结果 =====
st.download_button(
    "📥 下载 LogP 预测结果 CSV",
    df_valid.to_csv(index=False).encode("utf-8"),
    file_name="logp_predictions_xgb.csv",
    mime="text/csv"
)