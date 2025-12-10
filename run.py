"""
Advanced Composite Property Calculator - Batch Run Script

다양한 케이스를 설정하여 한번에 실행하고 결과를 CSV로 저장
"""

import csv
import time
from datetime import datetime
from effective_property_advanced import AdvancedCompositeCalculator


# ============================================================================
# 설정
# ============================================================================

# 재료 물성 (고정)
MATRIX_PROPS = {'E': 3500, 'nu': 0.35, 'alpha': 60e-6}
FIBER_PROPS = {'E': 72000, 'nu': 0.22, 'alpha': 5e-6}

# 해석 파일(.rst) 저장 여부
SAVE_RESULT_FILES = False

# 해석 케이스 리스트: 각 케이스는 딕셔너리로 정의
# 모든 키는 선택 키 (기본값: fiber_vf=0.1, ele_size=0.1, element_type='SOLID185')
# 선택 키: rve_size (기본값 1.0), strain_mag (기본값 0.001)
CASES = [
    # 기본 케이스 (모든 기본값 사용)
    {},

    # ele_size 변화
    {'ele_size': 0.04},
    {'ele_size': 0.06},
    {'ele_size': 0.08},
    {'ele_size': 0.12},

    # element_type 비교
    {'element_type': 'SOLID187'},

    # fiber_vf 변화
    {'fiber_vf': 0.20},
    {'fiber_vf': 0.30},
]

# 출력 CSV 파일명
OUTPUT_CSV = f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"


# ============================================================================
# 메인 실행
# ============================================================================

def run_all_cases():
    """모든 케이스를 순차적으로 실행하고 결과를 CSV로 저장"""

    results = []
    total_start = time.time()

    print(f"\n{'='*70}")
    print(f"BATCH ANALYSIS: {len(CASES)} cases")
    print(f"Save result files: {SAVE_RESULT_FILES}")
    print(f"{'='*70}")

    for i, case in enumerate(CASES):
        # 기본값 적용
        fiber_vf = case.get('fiber_vf', 0.1)
        ele_size = case.get('ele_size', 0.1)
        element_type = case.get('element_type', 'SOLID185')
        rve_size = case.get('rve_size', 1.0)
        strain_mag = case.get('strain_mag', 0.001)

        print(f"\n{'#'*70}")
        print(f"CASE {i+1}/{len(CASES)}: Vf={fiber_vf}, ele_size={ele_size}, {element_type}")
        print(f"{'#'*70}")

        case_start = time.time()

        # 계산기 인스턴스 생성
        calc = AdvancedCompositeCalculator(rve_size=rve_size, fiber_vf=fiber_vf)
        calc.set_material('matrix', **MATRIX_PROPS)
        calc.set_material('fiber', **FIBER_PROPS)

        try:
            # MAPDL 실행
            calc.launch(
                nproc=12,
                run_location='/tmp/mapdl_adv',
                override=True,
                additional_switches='-smp'
            )

            # 전체 해석 실행
            props = calc.run_full_analysis(
                ele_size=ele_size,
                strain_mag=strain_mag,
                element_type=element_type
            )

            # 결과 출력
            calc.print_results()

            # 해석 시간 계산
            elapsed = time.time() - case_start

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
                'nu_yz': props.get('nu_yz', ''),
                'nu_zx': props.get('nu_zx', ''),
                'alpha_x': props.get('alpha_x', ''),
                'alpha_y': props.get('alpha_y', ''),
                'alpha_z': props.get('alpha_z', ''),
                'time_sec': round(elapsed, 1),
                'status': 'OK'
            }
            results.append(result)

            print(f"\n  Case {i+1} completed in {elapsed:.1f} seconds")

        except Exception as e:
            elapsed = time.time() - case_start
            print(f"ERROR in case {i+1}: {e}")
            result = {
                'case': i + 1,
                'fiber_vf': fiber_vf,
                'rve_size': rve_size,
                'ele_size': ele_size,
                'element_type': element_type,
                'time_sec': round(elapsed, 1),
                'status': f'ERROR: {e}'
            }
            results.append(result)

        finally:
            # 해석 파일 삭제 (저장하지 않는 경우)
            if not SAVE_RESULT_FILES and calc.mapdl:
                cleanup_result_files(calc.mapdl.directory)
            calc.exit()

    # CSV 저장
    save_results_to_csv(results)

    total_elapsed = time.time() - total_start
    print(f"\n{'='*70}")
    print(f"BATCH COMPLETE: {len(results)} cases processed")
    print(f"Total time: {total_elapsed:.1f} seconds")
    print(f"Results saved to: {OUTPUT_CSV}")
    print(f"{'='*70}")

    return results


def cleanup_result_files(work_dir):
    """해석 결과 파일(.rst) 삭제"""
    import os
    import glob
    for filepath in glob.glob(os.path.join(work_dir, '*.rst')):
        try:
            os.remove(filepath)
        except OSError:
            pass


def save_results_to_csv(results):
    """결과를 CSV 파일로 저장"""
    if not results:
        return

    fieldnames = [
        'case', 'fiber_vf', 'rve_size', 'ele_size', 'element_type',
        'Ex', 'Ey', 'Ez', 'Gxy', 'Gyz', 'Gxz',
        'nu_xy', 'nu_yz', 'nu_zx',
        'alpha_x', 'alpha_y', 'alpha_z',
        'time_sec', 'status'
    ]

    with open(OUTPUT_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = {k: result.get(k, '') for k in fieldnames}
            writer.writerow(row)

    print(f"\nResults saved to {OUTPUT_CSV}")


if __name__ == "__main__":
    run_all_cases()
