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

## 사용법

### 기본 실행

```python
from effective_property_calculator import CompositeEffectivePropertyCalculator

# 계산기 인스턴스 생성
calc = CompositeEffectivePropertyCalculator(rve_size=1.0, fiber_vf=0.10)

# 재료 물성 설정 (선택사항)
calc.set_matrix_properties(E=3500, nu=0.35, alpha=60e-6)  # Epoxy
calc.set_fiber_properties(E=230000, nu=0.20, alpha=-0.5e-6)  # Carbon fiber

# MAPDL 시작
calc.start_mapdl()

# 유효 물성 계산
props = calc.calculate_effective_properties(element_size=0.05)

# 결과 출력
calc.print_results()

# MAPDL 종료
calc.stop_mapdl()
```

### 명령줄 실행

```bash
python effective_property_calculator.py
```

## 방법론

### 균질화 이론

RVE에 6가지 변형 상태를 적용하여 유효 강성 행렬을 구합니다:

1. **단축 변형 εxx**: Ex, νxy, νxz 계산
2. **단축 변형 εyy**: Ey, νyx, νyz 계산
3. **단축 변형 εzz**: Ez, νzx, νzy 계산
4. **전단 변형 γxy**: Gxy 계산
5. **전단 변형 γyz**: Gyz 계산
6. **전단 변형 γzx**: Gzx 계산

열팽창계수는 균일 온도 변화(ΔT = 1°C)를 적용하여 각 방향의 변형률로부터 계산합니다.

### 경계조건

주기적 경계조건(Periodic Boundary Conditions)을 적용하여 무한 반복 구조를 모사합니다.

## 기본 재료 물성

### Matrix (Epoxy)
- E = 3,500 MPa
- ν = 0.35
- α = 60 × 10⁻⁶ /°C

### Fiber (Carbon)
- E = 230,000 MPa
- ν = 0.20
- α = -0.5 × 10⁻⁶ /°C

## 출력 예시

```
============================================================
EFFECTIVE MATERIAL PROPERTIES
============================================================

--- Elastic Moduli (MPa) ---
  Ex = 5234.56
  Ey = 5234.56
  Ez = 26123.45

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
  CTEz = 1.23e-05
```

## 라이선스

MIT License
