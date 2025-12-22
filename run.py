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

# Reference 값 (직접 입력)
REFERENCE = {
    'Ex': 0.0,
    'Ey': 0.0,
    'Ez': 0.0,
    'Gxy': 0.0,
    'Gyz': 0.0,
    'Gxz': 0.0,
    'nu_xy': 0.0,
    'nu_yz': 0.0,
    'nu_xz': 0.0,
    'alpha_x': 0.0,
    'alpha_y': 0.0,
    'alpha_z': 0.0,
}

# 해석 케이스 리스트: 각 케이스는 딕셔너리로 정의
# 모든 키는 선택 키 (기본값: fiber_vf=0.1, ele_size=0.1, element_type='SOLID185')
# 선택 키: rve_size (기본값 1.0), strain_mag (기본값 0.001)
# DB 관련 키:
#   - save_db: DB 저장 파일명 (모델 빌드 후 저장)
#   - db_path: DB 불러오기 경로 (이 옵션 사용 시 모델 빌드 생략)
CASES = [
    # 기본 케이스 - 모델 빌드 후 DB 저장
    {'save_db': 'rve_vf10_ele01'},

    # DB에서 불러와서 해석 (모델 빌드 없이 바로 해석)
    # {'db_path': '/tmp/mapdl_adv/rve_vf10_ele01.db'},

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
        db_path = case.get('db_path', None)  # DB 불러오기 경로
        save_db = case.get('save_db', None)  # DB 저장 파일명

        # 케이스 정보 출력
        print(f"\n{'#'*70}")
        if db_path:
            print(f"CASE {i+1}/{len(CASES)}: [DB RESUME] {db_path}")
        else:
            print(f"CASE {i+1}/{len(CASES)}: Vf={fiber_vf}, ele_size={ele_size}, {element_type}")
            if save_db:
                print(f"  -> DB will be saved as: {save_db}.db")
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

            # DB resume 또는 새로 빌드
            result_prefix = f"case{i+1}_" if SAVE_RESULT_FILES else ''
            if db_path:
                # DB에서 불러와서 해석
                props = calc.run_analysis_from_db(
                    db_path=db_path,
                    strain_mag=strain_mag,
                    save_results=SAVE_RESULT_FILES,
                    result_prefix=result_prefix
                )
            else:
                # 새로 모델 빌드 후 해석
                props = calc.run_full_analysis(
                    ele_size=ele_size,
                    strain_mag=strain_mag,
                    element_type=element_type,
                    save_results=SAVE_RESULT_FILES,
                    result_prefix=result_prefix,
                    save_db=save_db
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
                'db_path': db_path if db_path else '',
                'Ex': props.get('Ex', ''),
                'Ey': props.get('Ey', ''),
                'Ez': props.get('Ez', ''),
                'Gxy': props.get('Gxy', ''),
                'Gyz': props.get('Gyz', ''),
                'Gxz': props.get('Gxz', ''),
                'nu_xy': props.get('nu_xy', ''),
                'nu_yz': props.get('nu_yz', ''),
                'nu_xz': props.get('nu_xz', ''),
                'alpha_x': props.get('alpha_x', ''),
                'alpha_y': props.get('alpha_y', ''),
                'alpha_z': props.get('alpha_z', ''),
                'time_sec': round(elapsed, 1),
                'status': 'OK'
            }

            # Reference와의 차이 계산
            for key in REFERENCE:
                ref_val = REFERENCE[key]
                calc_val = result.get(key, '')
                if ref_val != 0 and calc_val != '':
                    diff_pct = (calc_val - ref_val) / ref_val * 100
                    result[f'{key}_diff%'] = round(diff_pct, 2)
                else:
                    result[f'{key}_diff%'] = ''

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
                'db_path': db_path if db_path else '',
                'time_sec': round(elapsed, 1),
                'status': f'ERROR: {e}'
            }
            results.append(result)

        finally:
            calc.exit()

    # CSV 저장
    save_results_to_csv(results)

    # Plot 생성
    plot_results(results)

    total_elapsed = time.time() - total_start
    print(f"\n{'='*70}")
    print(f"BATCH COMPLETE: {len(results)} cases processed")
    print(f"Total time: {total_elapsed:.1f} seconds")
    print(f"Results saved to: {OUTPUT_CSV}")
    print(f"{'='*70}")

    return results


def save_results_to_csv(results):
    """결과를 CSV 파일로 저장"""
    if not results:
        return

    fieldnames = [
        'case', 'fiber_vf', 'rve_size', 'ele_size', 'element_type', 'db_path',
        'Ex', 'Ey', 'Ez', 'Gxy', 'Gyz', 'Gxz',
        'nu_xy', 'nu_yz', 'nu_xz',
        'alpha_x', 'alpha_y', 'alpha_z',
        'Ex_diff%', 'Ey_diff%', 'Ez_diff%',
        'Gxy_diff%', 'Gyz_diff%', 'Gxz_diff%',
        'nu_xy_diff%', 'nu_yz_diff%', 'nu_xz_diff%',
        'alpha_x_diff%', 'alpha_y_diff%', 'alpha_z_diff%',
        'time_sec', 'status'
    ]

    with open(OUTPUT_CSV, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            row = {k: result.get(k, '') for k in fieldnames}
            writer.writerow(row)

    print(f"\nResults saved to {OUTPUT_CSV}")


def plot_results(results):
    """결과 plot 생성"""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed. Skipping plot.")
        return

    # 성공한 케이스만 필터링
    ok_results = [r for r in results if r.get('status') == 'OK']
    if not ok_results:
        return

    # diff% 값이 있는 property만 plot
    diff_keys = ['Ex_diff%', 'Ey_diff%', 'Ez_diff%', 'Gxy_diff%', 'Gyz_diff%', 'Gxz_diff%',
                 'nu_xy_diff%', 'nu_yz_diff%', 'nu_xz_diff%',
                 'alpha_x_diff%', 'alpha_y_diff%', 'alpha_z_diff%']

    # 유효한 diff 데이터가 있는지 확인
    has_diff_data = any(r.get(k, '') != '' for r in ok_results for k in diff_keys)
    if not has_diff_data:
        print("No reference values set. Skipping diff plot.")
        return

    cases = [r['case'] for r in ok_results]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Results vs Reference (%diff)')

    # E moduli
    ax = axes[0, 0]
    for key in ['Ex_diff%', 'Ey_diff%', 'Ez_diff%']:
        vals = [r.get(key, None) for r in ok_results]
        if any(v is not None and v != '' for v in vals):
            ax.plot(cases, vals, 'o-', label=key.replace('_diff%', ''))
    ax.set_xlabel('Case')
    ax.set_ylabel('Diff (%)')
    ax.set_title('Elastic Moduli E')
    ax.legend()
    ax.grid(True)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

    # G moduli
    ax = axes[0, 1]
    for key in ['Gxy_diff%', 'Gyz_diff%', 'Gxz_diff%']:
        vals = [r.get(key, None) for r in ok_results]
        if any(v is not None and v != '' for v in vals):
            ax.plot(cases, vals, 'o-', label=key.replace('_diff%', ''))
    ax.set_xlabel('Case')
    ax.set_ylabel('Diff (%)')
    ax.set_title('Shear Moduli G')
    ax.legend()
    ax.grid(True)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

    # Poisson's ratio
    ax = axes[1, 0]
    for key in ['nu_xy_diff%', 'nu_yz_diff%', 'nu_xz_diff%']:
        vals = [r.get(key, None) for r in ok_results]
        if any(v is not None and v != '' for v in vals):
            ax.plot(cases, vals, 'o-', label=key.replace('_diff%', ''))
    ax.set_xlabel('Case')
    ax.set_ylabel('Diff (%)')
    ax.set_title("Poisson's Ratio")
    ax.legend()
    ax.grid(True)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

    # CTE
    ax = axes[1, 1]
    for key in ['alpha_x_diff%', 'alpha_y_diff%', 'alpha_z_diff%']:
        vals = [r.get(key, None) for r in ok_results]
        if any(v is not None and v != '' for v in vals):
            ax.plot(cases, vals, 'o-', label=key.replace('_diff%', ''))
    ax.set_xlabel('Case')
    ax.set_ylabel('Diff (%)')
    ax.set_title('Thermal Expansion Coeff.')
    ax.legend()
    ax.grid(True)
    ax.axhline(y=0, color='k', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    plot_file = OUTPUT_CSV.replace('.csv', '_plot.png')
    plt.savefig(plot_file, dpi=150)
    plt.close()
    print(f"Plot saved to {plot_file}")


if __name__ == "__main__":
    run_all_cases()
