import streamlit as st

st.set_page_config(page_title="2025CADD课程实践", page_icon="🧪", layout="wide")

# ==== 顶部标题 ====
st.markdown("""
    <div style="text-align: center;">
        <h1 style="color: #1F4E79; font-size: 40px;">🧬 2025 CADD课程实践平台</h1>
        <p style="color: #666; font-size: 18px;">
            一站式药物分子性质预测与评估系统，点击左侧模块进入各功能页面
        </p>
    </div>
""", unsafe_allow_html=True)

# ==== 样式定义 ====
st.markdown("""
    <style>
        .card {
            background: linear-gradient(145deg, #f4faff, #eaf4ff);
            border: 1px solid #d0e3f5;
            border-radius: 12px;
            padding: 20px;
            text-align: center;
            transition: all 0.3s ease;
            height: 150px;
        }
        .card:hover {
            transform: scale(1.02);
            box-shadow: 0 4px 12px rgba(0, 90, 180, 0.15);
        }
        .card-title {
            font-size: 20px;
            font-weight: bold;
            color: #1F4E79;
        }
        .card-desc {
            font-size: 14px;
            color: #333;
            margin-top: 8px;
        }
    </style>
""", unsafe_allow_html=True)

# ==== 一级模块：药物性评估入口 ====
st.subheader("🧪 药物性评估入口")
col1 = st.columns(1)[0]
with col1:
    st.markdown("""
        <div class="card">
            <div class="card-title">💊 药物性评估</div>
            <div class="card-desc">基于 QED 和 Lipinski 规则进行分子初筛，决定是否进入进一步性质预测。</div>
        </div>
    """, unsafe_allow_html=True)

# ==== 一级模块：三大核心性质预测 ====
st.subheader("🧬 分子性质预测模块")
col2, col3, col4 = st.columns(3)

with col2:
    st.markdown("""
        <div class="card">
            <div class="card-title">💧 溶解度预测</div>
            <div class="card-desc">使用机器学习预测 Log(Solubility)，评估分子水中溶解能力。</div>
        </div>
    """, unsafe_allow_html=True)
with col3:
    st.markdown("""
        <div class="card">
            <div class="card-title">🔥 脂溶性预测</div>
            <div class="card-desc">估算 LogP 值，反映分子穿膜能力及脂水分布倾向。</div>
        </div>
    """, unsafe_allow_html=True)
with col4:
    st.markdown("""
        <div class="card">
            <div class="card-title">🌊 水合自由能预测</div>
            <div class="card-desc">预测 ΔG_hydration 值，衡量分子与水相互作用强度。</div>
        </div>
    """, unsafe_allow_html=True)

# ==== 一级模块：知识拓展 ====
st.subheader("📚 拓展与参考模块")
col5 = st.columns(1)[0]
with col5:
    st.markdown("""
        <div class="card">
            <div class="card-title">📚 知识获取</div>
            <div class="card-desc">快速了解药物毒副作用、性质关系与生物活性机制。</div>
        </div>
    """, unsafe_allow_html=True)

# ==== 页脚 ====
st.markdown("""
    <hr>
    <div style="text-align: center; color: #aaa; font-size: 13px; margin-top: 30px;">
        © 2251110 宁书芮 计算机辅助药物设计课程实践与展示
    </div>
""", unsafe_allow_html=True)
