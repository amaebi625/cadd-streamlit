# train_hydration_model.py

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import joblib

# === 1. 加载数据 === #
df = pd.read_csv("data/SAMPL.csv")  # 请确保列名为 'smiles', 'expt'

# === 2. SMILES 转 Fingerprint === #
def smiles_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    else:
        return None

df["fp"] = df["smiles"].apply(smiles_to_fp)
df = df[df["fp"].notnull()].copy()
X = np.array(df["fp"].tolist())
y = df["expt"].values

# === 3. 拆分训练集/测试集 === #
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# === 4. 模型训练 === #
model = RandomForestRegressor(n_estimators=200, random_state=42)
model.fit(X_train, y_train)

# === 5. 性能评估 === #
y_pred = model.predict(X_test)
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)

r2 = r2_score(y_test, y_pred)

print(f"RMSE: {rmse:.3f} | R²: {r2:.3f}")

# === 6. 保存模型 === #
joblib.dump(model, "config/hydration_rf_model.pkl")
