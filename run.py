"""
Advanced Composite Property Calculator - Batch Run Script

다양한 케이스를 설정하여 한번에 실행하고 결과를 CSV로 저장
"""

import csv
from datetime import datetime
from effective_property_advanced import AdvancedCompositeCalculator


# ============================================================================
# 케이스 설정 (원하는 대로 수정)
# ============================================================================

# 재료 물성 (고정)
MATRIX_PROPS = {'E': 3500, 'nu': 0.35, 'alpha': 60e-6}
FIBER_PROPS = {'E': 72000, 'nu': 0.22, 'alpha': 5e-6}

# 해석 케이스 리스트: 각 케이스는 딕셔너리로 정의
# 필수 키: fiber_vf, ele_size, element_type
# 선택 키: rve_size (기본값 1.0), strain_mag (기본값 0.001)
CASES = [
    # fiber_vf = 0.1, 다양한 ele_size
    {'fiber_vf': 0.10, 'ele_size': 0.04, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.06, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.08, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.10, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.12, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.14, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.10, 'ele_size': 0.16, 'element_type': 'SOLID185'},

    # fiber_vf = 0.1, SOLID187 비교
    {'fiber_vf': 0.10, 'ele_size': 0.08, 'element_type': 'SOLID187'},

    # 다양한 fiber_vf
    {'fiber_vf': 0.20, 'ele_size': 0.08, 'element_type': 'SOLID185'},
    {'fiber_vf': 0.30, 'ele_size': 0.08, 'element_type': 'SOLID185'},
]

# 출력 CSV 파일명
OUTPUT_CSV = f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"


# ============================================================================
# 메인 실행
# ============================================================================

def run_all_cases():
    """모든 케이스를 순차적으로 실행하고 결과를 CSV로 저장"""

    results = []

    print(f"\n{'='*70}")
    print(f"BATCH ANALYSIS: {len(CASES)} cases")
    print(f"{'='*70}")

    for i, case in enumerate(CASES):
        fiber_vf = case['fiber_vf']
        ele_size = case['ele_size']
        element_type = case['element_type']
        rve_size = case.get('rve_size', 1.0)
        strain_mag = case.get('strain_mag', 0.001)

        print(f"\n{'#'*70}")
        print(f"CASE {i+1}/{len(CASES)}: Vf={fiber_vf}, ele_size={ele_size}, {element_type}")
        print(f"{'#'*70}")

        # 계산기 인스턴스 생성
        calc = AdvancedCompositeCalculator(rve_size=rve_size, fiber_vf=fiber_vf)
        calc.set_material('matrix', **MATRIX_PROPS)
        calc.set_material('fiber', **FIBER_PROPS)

        try:
            # MAPDL 실행
            calc.launch(run_location='/tmp/mapdl_adv', override=True, additional_switches='-smp')

            # 전체 해석 실행
            props = calc.run_full_analysis(
                ele_size=ele_size,
                strain_mag=strain_mag,
                element_type=element_type
            )

            # 결과 출력
            calc.print_results()

            # 결과 저장
            result = {
                'case': i + 1,
                'fiber_vf': fiber_vf,
                'rve_size': rve_size,
                'ele_size': ele_size,
                'element_type': element_type,
                'Ex': props.get('Ex', ''),
                'Ey': props.get('Ey', ''),
                'Ez': props.get('Ez', ''),
                'Gxy': props.get('Gxy', ''),
                'Gyz': props.get('Gyz', ''),
                'Gxz': props.get('Gxz', ''),
                'nu_xy': props.get('nu_xy', ''),
                'nu_yx': props.get('nu_yx', ''),
                'nu_xz': props.get('nu_xz', ''),
                'nu_zx': props.get('nu_zx', ''),
                'nu_yz': props.get('nu_yz', ''),
                'nu_zy': props.get('nu_zy', ''),
                'alpha_x': props.get('alpha_x', ''),
                'alpha_y': props.get('alpha_y', ''),
                'alpha_z': props.get('alpha_z', ''),
                'status': 'OK'
            }
            results.append(result)

        except Exception as e:
            print(f"ERROR in case {i+1}: {e}")
            result = {
                'case': i + 1,
                'fiber_vf': fiber_vf,
                'rve_size': rve_size,
                'ele_size': ele_size,
                'element_type': element_type,
                'status': f'ERROR: {e}'
            }
            results.append(result)

        finally:
            calc.exit()

    # CSV 저장
    save_results_to_csv(results)

    print(f"\n{'='*70}")
    print(f"BATCH COMPLETE: {len(results)} cases processed")
    print(f"Results saved to: {OUTPUT_CSV}")
    print(f"{'='*70}")

    return results


def save_results_to_csv(results):
    """결과를 CSV 파일로 저장"""
    if not results:
        return

    fieldnames = [
        'case', 'fiber_vf', 'rve_size', 'ele_size', 'element_type',
        'Ex', 'Ey', 'Ez', 'Gxy', 'Gyz', 'Gxz',
        'nu_xy', 'nu_yx', 'nu_xz', 'nu_zx', 'nu_yz', 'nu_zy',
        'alpha_x', 'alpha_y', 'alpha_z', 'status'
    ]

    with open(OUTPUT_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            # 누락된 필드는 빈 문자열로 채움
            row = {k: result.get(k, '') for k in fieldnames}
            writer.writerow(row)

    print(f"\nResults saved to {OUTPUT_CSV}")


if __name__ == "__main__":
    run_all_cases()
