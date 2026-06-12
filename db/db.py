ifrom sqlalchemy import create_engine
import streamlit as st

@st.cache_resource
def get_engine():
    return create_engine(
        st.secrets["DATABASE_URL"]
    )
