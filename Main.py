import streamlit as st

st.set_page_config(
    page_title="Cancer Methylation Explorer",
    page_icon="🧬",
    layout="wide"
)

st.title("Cancer Methylation Explorer")

st.markdown("""
Welcome to the methylation exploration portal.

Use the side menu to navigate between modules:

- Region Explorer
- Gene Explorer

""")
