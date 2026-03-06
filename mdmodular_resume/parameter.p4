# ===== P4: Continue / Resume / Slice Mode =====
# ZH: 该文件控制是新算、续算，还是从指定时间切片后执行注入-only模式。
# EN: This file controls fresh run, normal resume, or slice-and-run injection-only mode.
# JA: このファイルは新規計算、通常再開、または時刻スライス後の注入専用モードを制御します。

# --- Base Resume Switch / 基础续算开关 / 基本再開スイッチ ---
# ZH: 0=从 First_Data_file 新算并重置初速度；1=从 Data/Et0 续算。
#      =1 时，程序会自动扫描 FileName.txt 中 Data_file/Et_file 模板对应的全部编号文件（00/01/02...）。
# EN: 0=start fresh from First_Data_file with reinitialized velocities; 1=resume from Data/Et0.
#      =1, the program auto-scans all indexed files (00/01/02...) matching Data_file/Et_file templates in FileName.txt.
# JA: 0=First_Data_file から新規計算して初速度を再生成; 1=Data/Et0 から再開。
#      =1 の場合、FileName.txt の Data_file/Et_file テンプレートに一致する全番号ファイル（00/01/02...）を自動走査します。
continue_run = 1

# ZH: 续算时 0=取最后一帧；1=按 resume_time_fs 寻找时间点。
# EN: For resume, 0=use last frame; 1=use resume_time_fs.
# JA: 再開時、0=最終フレーム使用; 1=resume_time_fs を使用。
resume_use_target_time = 1

# ZH: 续算目标时间 (fs)。
# EN: Resume target time (fs).
# JA: 再開目標時刻 (fs)。
#resume_time_fs = -1.0
resume_time_fs = 70800

# ZH: 1=要求 Data 与 Et0 在目标时刻精确匹配（容差 0.01*dt）；0=允许最近点。
# EN: 1=require exact Data/Et0 match at target time (tol 0.01*dt); 0=allow nearest.
# JA: 1=目標時刻で Data/Et0 の厳密一致（許容 0.01*dt）; 0=最近点を許可。
resume_require_exact_time = 1

# ZH: 续算输出 1=自动写到下一个编号，0=覆盖当前 Data_file/Et_file。
# EN: Resume output 1=write to next index, 0=overwrite current Data_file/Et_file.
# JA: 再開出力 1=次番号へ書き込み, 0=現在の Data_file/Et_file を上書き。
continue_write_next_index = 1
###########################################################################################
# --- Resume Action / 续算动作 / 再開アクション ---
# ZH: 0=普通续算；1=切片到新 DataN 目录并进入注入-only模式。
# EN: 0=normal continue; 1=slice into new DataN folder and run injection-only mode.
# JA: 0=通常再開; 1=新しい DataN フォルダへスライスして注入専用モード実行。
resume_action = 1

# ZH: 当 resume_action=1 时，1=运行时询问切片时间；0=使用下方 slice_time_fs。
# EN: When resume_action=1, 1=ask slice time interactively; 0=use slice_time_fs below.
# JA: resume_action=1 のとき、1=実行時にスライス時刻を入力; 0=下の slice_time_fs を使用。
prompt_slice_time = 1

# ZH: 切片目标时间 (fs)，<0 表示最后一帧。
# EN: Slice target time (fs), <0 means last frame.
# JA: スライス目標時刻 (fs)、<0 は最終フレーム。
slice_time_fs = -1.0

# ZH: 切片模式下程序会自动把切片帧时间重置为 0 再开始计算（固定行为）。
# EN: In slice mode, the selected frame time is always reset to 0 before running (fixed behavior).
# JA: スライスモードでは、選択フレーム時刻は常に 0 にリセットしてから計算を開始します（固定動作）。

# ZH: 1=取最近时间切片；0=要求精确匹配。
# EN: 1=use nearest-time slice; 0=require exact match.
# JA: 1=最近時刻でスライス; 0=厳密一致を要求。
slice_use_nearest_time = 1

# ZH: 切片输出目录前缀，会自动生成 Data1/Data2/...。
# EN: Output directory prefix for slice mode, auto-creates Data1/Data2/....
# JA: スライスモードの出力ディレクトリ接頭辞。Data1/Data2/... を自動生成。
slice_output_prefix = Data

# ZH: 注入-only模式的时长与控温参数已移动到 parameter.p5。
# EN: Inject-only duration/thermostat parameters were moved to parameter.p5.
# JA: 注入専用モードの継続時間・温調パラメータは parameter.p5 へ移動しました。

# --- Force Model / 力模型 / 力モデル ---
# ZH: 当前支持 BeW_ABOP_LCL。
# EN: Currently supported: BeW_ABOP_LCL.
# JA: 現在対応: BeW_ABOP_LCL。
force_model = BeW_ABOP_LCL
