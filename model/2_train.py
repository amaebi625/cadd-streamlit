# 文件名：train_logp_model.py
import pandas as pd
import numpy as np
import joblib
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from xgboost import XGBRegressor

# 读取数据
data = pd.read_csv("data/Lipophilicity.csv")  # 包含 'smiles' 和 'exp' 列

def mol_to_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    arr = np.zeros((2048,))
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

# 分子指纹转化
fps = []
exps = []
for i, row in data.iterrows():
    fp = mol_to_fp(row['smiles'])
    if fp is not None:
        fps.append(fp)
        exps.append(row['exp'])

X = np.array(fps)
y = np.array(exps)

# 划分数据
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# 使用 XGBoost 模型
model = XGBRegressor(
    n_estimators=300,
    max_depth=6,
    learning_rate=0.05,
    subsample=0.8,
    colsample_bytree=0.8,
    random_state=42,
    n_jobs=-1,
    verbosity=1
)
model.fit(X_train, y_train)

# 模型评估
y_pred = model.predict(X_test)
print("R² Score:", r2_score(y_test, y_pred))
mse = mean_squared_error(y_test, y_pred)
rmse = np.sqrt(mse)
print("RMSE:", rmse)

# 保存模型
joblib.dump(model, "config/logp_xgb_model.pkl")
print("✅ 模型已保存为 config/logp_xgb_model.pkl")
