import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

# 웹페이지 기본 설정
st.set_page_config(page_title="토질역학 토압 계산기", layout="wide")

st.title("🧱 옹벽 주동토압 즉시 계산 및 시각화 프로그램")
st.write("교수님 돌발 문제 풀이 및 검산용 통합 시스템 (Rankine 이론 적용)")

# 사이드바: 입력창
st.sidebar.header("📋 문제 조건 입력")
H = st.sidebar.number_input("1. 옹벽 총 높이 H (m)", min_value=1.0, max_value=30.0, value=6.0, step=0.5)
gamma = st.sidebar.number_input("2. 흙 단위중량 gamma (kN/m³)", min_value=10.0, max_value=25.0, value=18.0, step=0.5)
phi = st.sidebar.number_input("3. 내부마찰각 phi (도)", min_value=0.0, max_value=50.0, value=30.0, step=1.0)
c = st.sidebar.number_input("4. 점착력 c (kN/m²)", min_value=0.0, max_value=100.0, value=0.0, step=1.0)
q = st.sidebar.number_input("5. 상재하중 q (kN/m²)", min_value=0.0, max_value=200.0, value=10.0, step=1.0)

water_yn = st.sidebar.radio("6. 지하수위 존재 여부", ["없음 (n)", "있음 (y)"])

if water_yn == "있음 (y)":
    zw = st.sidebar.number_input("    -> 지표면으로부터 지하수위 깊이 Zw (m)", min_value=0.0, max_value=H, value=3.0, step=0.5)
    gamma_sat = st.sidebar.number_input("    -> 흙의 포화단위중량 gamma_sat (kN/m³)", min_value=10.0, max_value=30.0, value=20.0, step=0.5)
else:
    zw = H
    gamma_sat = gamma

# --- 연산 루틴 (인장균열 및 교수님 돌발문제 완벽 검산용 버전) ---
phi_rad = math.radians(phi)
Ka = (1 - math.sin(phi_rad)) / (1 + math.sin(phi_rad))
gamma_w = 9.81

# 1. 인장균열 깊이 (zc) 계산
if c > 0 and Ka > 0:
    zc_calc = (2 * c - q * math.sqrt(Ka)) / (gamma * math.sqrt(Ka))
    zc = max(0.0, min(zc_calc, H))
else:
    zc = 0.0

# 2. 지점별 '실제' 토압 연산 (인장력을 살려두기 위해 max 제거)
p_tot_0 = q * Ka - 2 * c * math.sqrt(Ka)
p_tot_zw = (q + gamma * zw) * Ka - 2 * c * math.sqrt(Ka)

gamma_sub = gamma_sat - gamma_w
sigma_v_eff_H = q + gamma * zw + gamma_sub * (H - zw)
p_soil_H = sigma_v_eff_H * Ka - 2 * c * math.sqrt(Ka)
u_H = gamma_w * (H - zw) if H > zw else 0.0
p_tot_H = p_soil_H + u_H

# 3. 합력 및 작용점 계산 (교수님 정답 기준: 인장균열 구간 무시하고 플러스 면적만 계산)
if zc < zw:
    p_start = max(0.0, p_tot_0) if zc == 0 else 0.0
    P1 = ((p_start + p_tot_zw) / 2.0) * (zw - zc)
    y1 = (H - zw) + ((zw - zc) / 3.0) * ((2 * p_start + p_tot_zw) / (p_start + p_tot_zw)) if (p_start + p_tot_zw) > 0 else (H - zw)
    
    P2 = ((max(0.0, p_tot_zw) + max(0.0, p_tot_H)) / 2.0) * (H - zw)
    y2 = ((H - zw) / 3.0) * ((2 * max(0.0, p_tot_zw) + max(0.0, p_tot_H)) / (max(0.0, p_tot_zw) + max(0.0, p_tot_H))) if (max(0.0, p_tot_zw) + max(0.0, p_tot_H)) > 0 else 0.0
else:
    P1 = 0.0
    y1 = 0.0
    P2 = ((0.0 + max(0.0, p_tot_H)) / 2.0) * (H - zc)
    y2 = ((H - zc) / 3.0)
    
Pa_total = P1 + P2
y_bar_total = (P1 * y1 + P2 * y2) / Pa_total if Pa_total > 0 else 0.0

# 4. 표와 그래프 출력을 위한 시각화용 보정 변수
p_tot_0_disp = max(0.0, p_tot_0)
p_tot_zw_disp = max(0.0, p_tot_zw)
p_tot_H_disp = max(0.0, p_tot_H)

# --- 레이아웃 배치 ---
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📊 연산 결과 및 검산 조견표")
    st.metric(label="주동토압계수 (Ka)", value=f"{Ka:.4f}")
    
    df_res = pd.DataFrame({
        "위치 (지표면 기준 깊이)": [f"z = 0m (지표면)", f"z = {zw}m (수위면)", f"z = {H}m (옹벽저면)"],
        "총 수평토압 (kN/m²)": [f"{p_tot_0_disp:.2f}", f"{p_tot_zw_disp:.2f}", f"{p_tot_H_disp:.2f}"],
        "간극수압 성분 (kN/m²)": ["0.00", "0.00", f"{u_H:.2f}"]
    })
    st.table(df_res)
    
    st.success(f"💥 **총 주동토압 합력 (Pa)** = **{Pa_total:.2f} kN/m**")
    st.info(f"📍 **토압 작용위치 (도심)** = 바닥에서 위로 **{y_bar_total:.2f} m** 지점")

with col2:
    st.subheader("📈 옹벽 배면 토압 분포도 그림")
    
    fig, ax = plt.subplots(figsize=(5, 6))
    z_layers = [0, zw, H]
    p_layers = [p_tot_0_disp, p_tot_zw_disp, p_tot_H_disp]
    
    ax.plot(p_layers, z_layers, color="red", linewidth=2, marker="o", label="Total Earth Pressure")
    ax.fill_betweenx(z_layers, 0, p_layers, color="orange", alpha=0.3)
    
    ax.set_xlim(left=0)
    ax.set_ylim(H, 0) # 깊이 방향 아래로 가도록 역축 설정
    ax.set_xlabel("Lateral Earth Pressure (kN/m²)")
    ax.set_ylabel("Depth z (m)")
    ax.grid(True, linestyle="--", alpha=0.6)
    
    # 작용점 위치 표시 화살표
    if Pa_total > 0:
        ax.annotate(f"Resultant Pa\n(y={y_bar_total:.2f}m)", xy=(Pa_total/H, H - y_bar_total), xytext=(Pa_total/H + 5, H - y_bar_total - 0.5),
                    arrowprops=dict(facecolor='blue', shrink=0.05, width=1.5, headwidth=6))
        
    st.pyplot(fig)