# Composite Effective Property Calculator

PyANSYS (PyMAPDL)를 이용한 복합재료 유효 물성 계산기

## 개요

이 프로그램은 Representative Volume Element (RVE) 균질화 방법을 사용하여 섬유 강화 복합재료의 유효 물성을 계산합니다.

### 형상
- **RVE 크기**: 1 x 1 x 1 mm 정육면체
- **Matrix**: 정육면체의 외곽 영역
- **Fiber**: 중심에 위치한 사각기둥 섬유 (Z축 방향)
- **섬유 체적분율**: 10%

### 계산 물성
- **탄성계수**: Ex, Ey, Ez (MPa)
- **전단계수**: Gxy, Gyz, Gzx (MPa)
- **푸아송비**: νxy, νyz, νzx
- **열팽창계수**: CTEx, CTEy, CTEz (1/°C)

## 요구사항

- Python 3.8+
- ANSYS MAPDL 라이선스
- PyMAPDL (`ansys-mapdl-core`)

## 설치

```bash
pip install -r requirements.txt
```

## 파일 구조

```
├── effective_property_calculator.py  # 기본 버전
├── effective_property_advanced.py    # 고급 버전 (Periodic BC + Hex mesh)
├── requirements.txt
└── README.md
```

## 실행 방법

### 방법 1: 기본 버전 (effective_property_calculator.py)

```bash
python effective_property_calculator.py
```

### 방법 2: 고급 버전 (effective_property_advanced.py) - 권장

고급 버전은 ANSYS Material Designer와 동일한 방식의 **Periodic Boundary Conditions**을 사용합니다.

```bash
python effective_property_advanced.py
```

또는 Python 코드에서:

```python
from effective_property_advanced import AdvancedCompositeCalculator

# 계산기 인스턴스 생성
calc = AdvancedCompositeCalculator(rve_size=1.0, fiber_vf=0.10)

# 재료 물성 설정
calc.set_material('matrix', E=3500, nu=0.35, alpha=60e-6)
calc.set_material('fiber', E=72000, nu=0.22, alpha=5e-6)

# MAPDL 시작
calc.launch(run_location='/tmp/mapdl_adv', override=True)

# 전체 해석 실행
props = calc.run_full_analysis(n_div=10, strain_mag=0.001)

# 결과 출력
calc.print_results()

# MAPDL 종료
calc.exit()
```

## 고급 버전의 주요 특징

### 1. 사각기둥 섬유 형상 (Square Fiber)
- 원통형 대신 사각기둥 형상 사용
- 대향면(opposite faces)의 메시가 완벽히 일치 가능
- Periodic boundary condition 적용에 필수

### 2. SOLID185 Hex 메시
- 8절점 육면체 요소(Hexahedral element) 사용
- Mapped meshing으로 정렬된 메시 생성
- 대향면 노드들이 1:1 대응

### 3. Periodic Boundary Conditions (CE 명령어)
ANSYS Material Designer와 동일한 방식의 주기적 경계조건 구현:

```
u(x⁺) - u(x⁻) = ε̄ · Δx
```

- **Master Node 방식**: 3개의 마스터 노드로 macroscopic strain 제어
- **CE (Constraint Equation)** 명령어로 대향면 노드 커플링
- 과구속(over-constraint) 없이 자연스러운 변형 허용

### 4. 경계조건 원리

각 대향면 쌍에 대해:
- X면 (X=L ↔ X=0): `u(L,y,z) - u(0,y,z) = [ε₁₁, ε₁₂, ε₃₁]ᵀ × L`
- Y면 (Y=L ↔ Y=0): `u(x,L,z) - u(x,0,z) = [ε₁₂, ε₂₂, ε₂₃]ᵀ × L`
- Z면 (Z=L ↔ Z=0): `u(x,y,L) - u(x,y,0) = [ε₃₁, ε₂₃, ε₃₃]ᵀ × L`

## RST 결과 파일

각 load case별로 ANSYS RST (result) 파일이 저장됩니다:

| Load Case | 파일명 | 설명 |
|-----------|--------|------|
| LC1 | `LC1_e11.rst` | 단축 변형 ε₁₁ |
| LC2 | `LC2_e22.rst` | 단축 변형 ε₂₂ |
| LC3 | `LC3_e33.rst` | 단축 변형 ε₃₃ |
| LC4 | `LC4_g12.rst` | 전단 변형 γ₁₂ |
| LC5 | `LC5_g23.rst` | 전단 변형 γ₂₃ |
| LC6 | `LC6_g31.rst` | 전단 변형 γ₃₁ |
| LC7 | `LC7_thermal.rst` | 열팽창 (ΔT=1°C) |

## 방법론

### 균질화 이론

RVE에 6가지 변형 상태를 적용하여 유효 강성 행렬을 구합니다:

1. **단축 변형 ε₁₁**: 6x6 강성행렬의 1열 계산
2. **단축 변형 ε₂₂**: 6x6 강성행렬의 2열 계산
3. **단축 변형 ε₃₃**: 6x6 강성행렬의 3열 계산
4. **전단 변형 γ₁₂**: 6x6 강성행렬의 4열 계산
5. **전단 변형 γ₂₃**: 6x6 강성행렬의 5열 계산
6. **전단 변형 γ₃₁**: 6x6 강성행렬의 6열 계산

강성행렬(C)의 역행렬인 유연성행렬(S)로부터 공학적 상수 추출:
- `S₁₁ = 1/Ex`, `S₂₂ = 1/Ey`, `S₃₃ = 1/Ez`
- `S₄₄ = 1/Gxy`, `S₅₅ = 1/Gyz`, `S₆₆ = 1/Gzx`
- `S₁₂ = -νxy/Ex`, etc.

### 열팽창계수

주기적 경계조건 하에서 균일 온도 변화(ΔT = 1°C)를 적용:
- 마스터 노드의 법선방향 변위만 자유롭게 두고 전단 변형 구속
- 마스터 노드 변위로부터 직접 열변형률 계산

## 기본 재료 물성

### Matrix (Epoxy)
- E = 3,500 MPa
- ν = 0.35
- α = 60 × 10⁻⁶ /°C

### Fiber (Glass)
- E = 72,000 MPa
- ν = 0.22
- α = 5 × 10⁻⁶ /°C

## 참고문헌

1. Xia, Z., Zhou, C., Yong, Q., Wang, X. (2006). "On selection of repeated unit cell model and application of unified periodic boundary conditions"
2. ANSYS Material Designer Theory Guide
3. Michael Okereke, Simeon Keates (2018). "Finite Element Applications: A Practical Guide to the FEM Process"

## 라이선스

MIT License
