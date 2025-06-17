import streamlit as st
import pandas as pd
import os, re
from io import StringIO
from Bio import Entrez
import xml.etree.ElementTree as ET

# ===== 设置 Entrez 用户邮箱 =====
Entrez.email = "your_email@example.com"

# ===== 页面标题与说明 =====
st.set_page_config(page_title="📚 知识获取", page_icon="📖", layout="wide")
st.title("📚 知识获取：毒副作用文献检索与信息提取")
st.markdown("通过 PubMed Central 检索临床毒理文献，并结合 GPT 提取相关毒副作用信息")

# ===== 搜索输入 =====
keyword = st.text_input("🔍 输入检索关键词", '"Clinical Toxicology" AND "Chemical"')

# ===== 执行搜索 =====
if keyword:
    st.markdown("正在搜索 PMC 文献...")
    try:
        search_handle = Entrez.esearch(db="pmc", term=keyword, retmode="xml", retmax=5)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        id_list = search_record['IdList']
    except Exception as e:
        st.error(f"❌ PubMed 搜索失败: {e}")
        id_list = []

    if id_list:
        st.success(f"共返回 {len(id_list)} 条文献（最多显示前1篇）")
        pmcid = id_list[0]

        try:
            fetch_handle = Entrez.efetch(db="pmc", id=pmcid, retmode="xml")
            xml_data = fetch_handle.read()
            fetch_handle.close()

            # 使用 ElementTree 解析全文
            root = ET.fromstring(xml_data)

            # ==== 提取标题 ====
            title = ""
            for tag in root.iter("article-title"):
                title = "".join(tag.itertext()).strip()
                break
            st.info(f"📖 文献标题：{title if title else '未找到标题'}")

            # ==== 提取摘要 ====
            abstract = ""
            for abstract_tag in root.iter("abstract"):
                parts = []
                for p in abstract_tag.iter("p"):
                    parts.append("".join(p.itertext()).strip())
                abstract += "\n".join(parts).strip()
            if abstract:
                st.success("✅ 成功提取摘要")
                st.text_area("📄 摘要内容", abstract, height=150)
            else:
                st.warning("⚠️ 未找到摘要内容")

            # ==== 提取正文 ====
            full_text = ""
            for sec in root.iter("sec"):
                for p in sec.iter("p"):
                    full_text += re.sub(r"<.*?>", "", "".join(p.itertext()).strip()) + "\n"
            st.text_area("📑 正文内容（截断）", full_text[:2000], height=250)

        except Exception as e:
            st.error(f"❌ 文献解析失败：{e}")
            abstract = ""


            # ==== 正文提取 ====
            full_text = ""
            try:
                body_sec = record.get("body", {}).get("sec", [])
                for sec in body_sec:
                    ps = sec.get("p", [])
                    if isinstance(ps, list):
                        for p in ps:
                            full_text += re.sub(r"<.*?>", "", str(p)) + "\n"
                    elif isinstance(ps, str):
                        full_text += re.sub(r"<.*?>", "", ps) + "\n"
                st.text_area("📑 正文内容（截断）", full_text[:2000], height=250)
            except:
                st.warning("⚠️ 正文内容提取失败")

            # ==== OpenAI 提取毒副作用信息 ====
            key = st.text_input("🔐 输入 OpenAI API Key", type="password")
            if key and abstract:
                try:
                    from openai import OpenAI
                    client = OpenAI(api_key=key)

                    # 普通 GPT 摘要
                    query = f"请从以下文献中提取与毒副作用相关的化合物信息，包括名字，类型和毒副作用描述：\n{abstract}"
                    response1 = client.chat.completions.create(
                        model="gpt-4",
                        messages=[{"role": "user", "content": query}]
                    )
                    gpt_output = response1.choices[0].message.content
                    st.text_area("🧠 GPT 提取结果（自然语言）", gpt_output, height=180)

                    # TSV 格式提取
                    tsv_prompt = f"""请从以下文献中提取与毒副作用相关的化合物信息，要求如下：
1. 仅输出获取的信息，不要输出额外的文字，英文回复；
2. 按照TSV格式输出，格式为："Compound\tType\tToxic Effect"；
3. 如果缺失信息请置空；
文献信息如下：
{abstract}"""

                    response2 = client.chat.completions.create(
                        model="gpt-4",
                        messages=[{"role": "user", "content": tsv_prompt}]
                    )
                    tsv_text = response2.choices[0].message.content.strip()
                    st.text_area("📄 GPT TSV 格式结果", tsv_text, height=160)

                    # 转为表格
                    try:
                        df = pd.read_csv(StringIO("Compound\tType\tToxic Effect\n" + tsv_text), sep="\t")
                        st.dataframe(df)
                    except:
                        st.warning("⚠️ GPT 返回格式无法正确解析为表格")

                except Exception as e:
                    st.error(f"❌ OpenAI API 失败：{e}")
