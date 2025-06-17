# train_esol_model.py
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error

# === 加载数据 ===
df = pd.read_csv("data/delaney-processed.csv")  # 请替换为你自己的数据文件名
df = df.dropna(subset=["smiles", "measured log solubility in mols per litre"])

# === 分子描述符函数 ===
def compute_features(smi):
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

features = df["smiles"].apply(compute_features)
features = features.dropna()
X = np.array(features.tolist())
y = df.loc[features.index, "measured log solubility in mols per litre"]

# === 划分数据，训练模型 ===
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# === 评估性能 ===
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
print(f"R² Score: {r2_score(y_test, y_pred):.4f}")
print(f"RMSE: ",rmse)

# === 保存模型 ===
joblib.dump(model, "config/esol_rf_model.pkl")
print("✅ 模型已保存到 config/esol_rf_model.pkl")
