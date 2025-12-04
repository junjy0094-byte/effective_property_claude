# Composite Effective Property Calculator

PyANSYS (PyMAPDL)를 이용한 복합재료 유효 물성 계산기

## 개요

이 프로그램은 Representative Volume Element (RVE) 균질화 방법을 사용하여 섬유 강화 복합재료의 유효 물성을 계산합니다.

### 형상
- **RVE 크기**: 1 x 1 x 1 mm 정육면체
- **Matrix**: 정육면체 전체
- **Fiber**: 중심에 위치한 원통형 섬유 (Z축 방향)
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
├── effective_property_advanced.py    # 고급 버전 (완전한 6x6 강성행렬 계산)
├── requirements.txt
└── README.md
```

## 실행 방법

### 방법 1: 기본 버전 (effective_property_calculator.py)

```bash
python effective_property_calculator.py
```

또는 Python 코드에서:

```python
from effective_property_calculator import CompositeEffectivePropertyCalculator

# 계산기 인스턴스 생성
calc = CompositeEffectivePropertyCalculator(rve_size=1.0, fiber_vf=0.10)

# 재료 물성 설정 (선택사항 - 기본값: Glass fiber + Epoxy)
calc.set_matrix_properties(E=3500, nu=0.35, alpha=60e-6)  # Epoxy
calc.set_fiber_properties(E=72000, nu=0.22, alpha=5e-6)   # Glass fiber

# MAPDL 시작
calc.start_mapdl(run_location='/tmp/mapdl_run', override=True)

# 유효 물성 계산
props = calc.calculate_effective_properties(element_size=0.05, strain_val=0.001)

# 결과 출력
calc.print_results()

# MAPDL 종료
calc.stop_mapdl()
```

### 방법 2: 고급 버전 (effective_property_advanced.py)

고급 버전은 완전한 6x6 강성행렬을 계산하고 역행렬(compliance matrix)로부터 공학적 상수를 추출합니다.

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
props = calc.run_full_analysis(elem_size=0.05, strain_mag=0.001)

# 결과 출력
calc.print_results()

# MAPDL 종료
calc.exit()
```

## RST 결과 파일

각 load case별로 ANSYS RST (result) 파일이 저장됩니다. 이 파일들은 ANSYS Mechanical 또는 다른 post-processor에서 열어 확인할 수 있습니다.

| Load Case | 파일명 | 설명 |
|-----------|--------|------|
| LC1 | `LC1_uniaxial_X.rst` / `LC1_e11.rst` | 단축 변형 εxx |
| LC2 | `LC2_uniaxial_Y.rst` / `LC2_e22.rst` | 단축 변형 εyy |
| LC3 | `LC3_uniaxial_Z.rst` / `LC3_e33.rst` | 단축 변형 εzz |
| LC4 | `LC4_shear_XY.rst` / `LC4_g12.rst` | 전단 변형 γxy |
| LC5 | `LC5_shear_YZ.rst` / `LC5_g23.rst` | 전단 변형 γyz |
| LC6 | `LC6_shear_ZX.rst` / `LC6_g31.rst` | 전단 변형 γzx |
| LC7 | `LC7_thermal.rst` | 열팽창 (ΔT=1°C) |

결과 파일 위치는 실행 시 출력됩니다:
```
Result files will be saved in: /tmp/mapdl_run
```

## 방법론

### 균질화 이론

RVE에 6가지 변형 상태를 적용하여 유효 강성 행렬을 구합니다:

1. **단축 변형 εxx**: C11, C21, C31 계산
2. **단축 변형 εyy**: C12, C22, C32 계산
3. **단축 변형 εzz**: C13, C23, C33 계산
4. **전단 변형 γxy**: C44 (Gxy) 계산
5. **전단 변형 γyz**: C55 (Gyz) 계산
6. **전단 변형 γzx**: C66 (Gzx) 계산

열팽창계수는 균일 온도 변화(ΔT = 1°C)를 적용하여 각 방향의 변형률로부터 계산합니다.

### 경계조건

Kinematic Uniform Boundary Conditions (KUBC)을 적용:
- 음의 면 (X-, Y-, Z-): 해당 방향 변위 고정
- 양의 면 (X+, Y+, Z+): 변형에 해당하는 변위 적용

## 기본 재료 물성

### Matrix (Epoxy)
- E = 3,500 MPa
- ν = 0.35
- α = 60 × 10⁻⁶ /°C

### Fiber (Glass)
- E = 72,000 MPa
- ν = 0.22
- α = 5 × 10⁻⁶ /°C

## 출력 예시

```
============================================================
COMPOSITE EFFECTIVE PROPERTY CALCULATION
============================================================
Result files will be saved in: /tmp/mapdl_run

Creating RVE geometry...
  RVE size: 1.0 x 1.0 x 1.0 mm
  Fiber radius: 0.1784 mm
  Fiber volume fraction: 10.0%

Load Case 1: Uniaxial strain εxx...
    Result saved: LC1_uniaxial_X.rst
...

============================================================
EFFECTIVE MATERIAL PROPERTIES
============================================================

--- Elastic Moduli (MPa) ---
  Ex = 4523.45
  Ey = 4523.45
  Ez = 10234.56

--- Shear Moduli (MPa) ---
  Gxy = 1456.78
  Gyz = 1823.45
  Gzx = 1823.45

--- Poisson's Ratios ---
  νxy = 0.3412
  νyz = 0.3156
  νzx = 0.3156

--- Coefficients of Thermal Expansion (1/°C) ---
  CTEx = 4.52e-05
  CTEy = 4.52e-05
  CTEz = 2.34e-05
```

## 라이선스

MIT License
