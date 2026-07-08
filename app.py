"""
CryoProtect-Screener: Web application for screening cryoprotectant proteins
for cryopreservation of donor organs and biomaterials.
"""

import streamlit as st
import requests
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import plotly.express as px
from io import StringIO

# ============================================================
# PAGE CONFIGURATION
# ============================================================
st.set_page_config(
    page_title="CryoProtect-Screener",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================
# LANGUAGE: All text strings for both languages
# ============================================================
TEXTS = {
    "ru": {
        "title": "🧬 КриоПротектор-Скринер",
        "subtitle": "*In silico* скрининг белков для криоконсервации биоматериала",
        "sidebar_header": "📋 О приложении",
        "sidebar_text": """**Проблема:** Каждые 17 минут умирает пациент, ожидающий трансплантации. Донорские органы хранятся всего 4–24 часа. Криоконсервация при -196°C могла бы продлить хранение до нескольких лет, создав «банки органов».

**Решение:** Антифризные белки (AFP), белки теплового шока (HSP) и LEA-белки — природные криопротекторы с низкой токсичностью.

**Как использовать:**
1. Введите UniProt ID белка
2. Нажмите «Анализировать»
3. Получите CryoScore и рекомендации по применению

**Автор:** Бондаренко Вероника""",
        "sidebar_footer": "Поддержите науку! Донорство органов спасает жизни.",
        "select_preset": "Или выберите из примеров:",
        "input_label": "UniProt ID:",
        "input_placeholder": "Например: P19479",
        "criteria_header": "🎯 Критерии криопротектора",
        "criteria_text": """**Универсальные критерии:**
1. ❄️ Связывает лёд/воду
2. 💧 Гидрофильный (GRAVY < 0)
3. 📏 Малый размер (< 30 кДа)
4. 🧬 Стабильный (индекс нестабильности < 40)
5. ✅ Безопасный (предсказанная нетоксичность)

**Специфические критерии:**
- HSP: термостабильность > 60°C
- LEA: высокая гидрофильность, защита от дегидратации
- Пептиды: способность проникать через клеточную мембрану""",
        "analyze_button": "🔬 Анализировать белок",
        "warning_id": "⚠️ Введите UniProt ID!",
        "loading_sequence": "Загрузка последовательности из UniProt...",
        "error_load": "❌ Не удалось загрузить",
        "check_id": "Проверьте ID.",
        "success_load": "✅ Белок загружен:",
        "type_hsp": "🔬 Тип белка: Белок теплового шока (HSP) — молекулярный шаперон, защита от стресса",
        "type_lea": "🔬 Тип белка: LEA-белок — защита от высыхания (ангидробиоз)",
        "type_afp": "🔬 Тип белка: Антифризный белок (AFP) — связывание льда, термический гистерезис",
        "type_other": "🔬 Тип белка: Синтетический/другой — требует экспериментальной валидации",
        "analyzing": "Анализ физико-химических свойств...",
        "error_short": "❌ Последовательность слишком короткая для анализа.",
        "results_header": "📊 Результаты анализа",
        "cryoscore_label": "CryoScore",
        "safety_label": "Безопасность",
        "penetration_label": "Проникновение в клетки",
        "excellent_delta": "Отличный кандидат!",
        "moderate_delta": "Средний",
        "low_delta": "Низкий приоритет",
        "passport_header": "📋 Паспорт кандидата в криопротекторы",
        "passport_param_col": "Параметр",
        "passport_value_col": "Значение",
        "passport_score_col": "Оценка",
        "passport_rows": ["UniProt ID", "Длина (а.о.)", "Мол. масса (кДа)", "GRAVY (гидрофобность)", "T+S+A% (лед-связывание)", "Изоэлектрическая точка", "Индекс нестабильности", "Токсичность", "Аллергенность", "Стабильность in vitro", "Проникновение в клетки"],
        "score_small": "✅ Малый", "score_large": "⚠️ Крупный",
        "score_light": "✅ Лёгкий", "score_heavy": "⚠️ Тяжёлый",
        "score_hydrophilic": "✅ Гидрофильный", "score_hydrophobic": "⚠️ Гидрофобный",
        "score_high": "✅ Высокий", "score_low": "⚠️ Низкий",
        "score_stable": "✅ Стабилен", "score_unstable": "⚠️ Нестабилен",
        "score_yes": "✅", "score_warn": "⚠️",
        "stable_text": "Стабилен", "unstable_text": "Нестабилен",
        "toxic_large": "Нетоксичен (крупный белок)",
        "toxic_short": "Возможна токсичность (короткий пептид)",
        "toxic_likely": "Вероятно, нетоксичен",
        "allergen_possible": "Возможна аллергенность",
        "allergen_likely": "Вероятно, неаллергенен",
        "safe_yes": "✅ БЕЗОПАСЕН",
        "safe_check": "⚠️ Требует проверки",
        "pen_strong": "✅ CPP+ (сильное)",
        "pen_moderate": "🟡 Умеренное",
        "pen_weak": "🟡 Слабое (малый размер)",
        "pen_none": "❌ Не проникает",
        "cryoscore_chart": "📈 Компоненты CryoScore",
        "comp_ice": "Лед-связывание\n(35%)",
        "comp_hydro": "Гидрофильность\n(25%)",
        "comp_size": "Размер\n(20%)",
        "comp_hbond": "H-связи\n(10%)",
        "comp_stab": "Стабильность\n(10%)",
        "chart_x": "Компонент", "chart_y": "Баллы",
        "chart_cryo_title": "Вклад компонентов в CryoScore",
        "kd_header": "💧 Профиль гидрофобности (Kyte-Doolittle)",
        "kd_x": "Позиция (окно 5 а.о.)", "kd_y": "Гидрофобность",
        "kd_title": "Скользящая гидрофобность (окно 5 а.о.)",
        "kd_boundary": "Граница гидрофобности",
        "aa_header": "🧬 Аминокислотный состав",
        "aa_thr": "Треонин (T)", "aa_ser": "Серин (S)", "aa_ala": "Аланин (A)",
        "aa_arg": "Аргинин (R)", "aa_lys": "Лизин (K)", "aa_other": "Прочие",
        "aa_title": "Состав белка",
        "func_hydrophobic": "Гидрофобные", "func_polar": "Полярные (T,S)",
        "func_charged": "Заряженные (R,K)", "func_other": "Прочие",
        "func_title": "Функциональные группы",
        "aa_col": "Аминокислота", "pct_col": "Процент",
        "group_col": "Группа",
        "rec_header": "💡 Рекомендация",
        "rec_excellent_title": "🏆 ОТЛИЧНЫЙ КАНДИДАТ!",
        "rec_moderate_title": "🟡 СРЕДНИЙ КАНДИДАТ",
        "rec_low_title": "🔴 НИЗКИЙ ПРИОРИТЕТ",
        "rec_moderate_text": "Можно тестировать, но не в первую очередь. Рекомендуется сравнить с другими белками из библиотеки.",
        "rec_low_text": "Рекомендуется поискать другие белки. Обратите внимание на белки с высоким T+S+A% и отрицательным GRAVY.",
        "rec_hsp": """**Рекомендации для криоконсервации спермы:**
1. Добавление в криопротекторную среду (0.1–1 мг/мл)
2. Инкубация со сперматозоидами (30 мин, 37°C)
3. Заморозка в парах жидкого азота
4. Оценка подвижности после разморозки (CASA)""",
        "rec_lea": """**Рекомендации для сухой консервации:**
1. Электропорация белка в клетки-мишени
2. Высушивание при 60% относительной влажности (24 ч)
3. Хранение при комнатной температуре
4. Регидратация и оценка выживаемости (MTT-тест)""",
        "rec_afp": """**Рекомендации для криоконсервации органов:**
1. Синтез рекомбинантного белка в *E. coli*
2. Измерение термического гистерезиса (наноосмометр)
3. MTT-тест на культуре HEK293
4. Криоконсервация гепатоцитов (0.1–1 мг/мл)
5. Сравнение выживаемости с DMSO-контролем""",
        "download_header": "📥 Скачать результаты",
        "download_button": "📄 Скачать паспорт (CSV)",
        "kb_header": "📚 База знаний",
        "kb_afp_title": "ℹ️ Что такое антифризные белки (AFP)?",
        "kb_afp_text": """**Антифризные белки (AFP)** — природные криопротекторы, обнаруженные у организмов, обитающих при отрицательных температурах (рыбы, насекомые, растения, бактерии).

**Механизм действия:**
- AFP адсорбируются на поверхности микроскопических кристаллов льда
- Подавляют их рост (эффект термического гистерезиса)
- Понижают точку замерзания, не изменяя точку плавления

**Классификация AFP:**
- **Type I** (рыбы): альфа-спиральные, обогащены аланином
- **Type II** (рыбы): C-type лектины, бета-сэндвич
- **Type III** (рыбы): глобулярные, бета-сэндвич
- **AFGP** (рыбы): гликопротеины с повторами Thr-Ala-Ala
- **Насекомые:** гиперактивные, правозакрученный бета-соленоид
- **Растительные:** хитиназы, тауматин-подобные белки, LTP
- **Бактериальные:** RTX-повторы, Ca²⁺-зависимые адгезины

**Эволюция:** AFP возникали конвергентно минимум 6–7 раз в разных таксонах.""",
        "kb_cryo_title": "🧮 Как рассчитывается CryoScore?",
        "kb_cryo_text": """**Формула:** CryoScore = 0.35×TSA + 0.25×(-GRAVY) + 0.20×(1000/MW) + 0.10×HBonds + 0.10×Stability

**Компоненты:**
- **T+S+A% (35%):** лед-связывающий потенциал — Thr, Ser, Ala формируют поверхность, комплементарную льду
- **-GRAVY (25%):** гидрофильность — чем гидрофильнее белок, тем лучше растворимость в криопротекторном растворе
- **1/MW (20%):** размер — низкомолекулярные белки эффективнее проникают в ткани
- **H-связи (10%):** плотность водородных связей коррелирует со стабильностью структуры
- **Стабильность (10%):** инвертированный индекс нестабильности (< 40 = стабилен in vitro)

**Интерпретация:**
- > 70: отличный кандидат, рекомендуется тестировать in vitro
- 50–70: средний кандидат
- < 50: низкий приоритет""",
        "kb_lab_title": "🔬 Как проверяют криопротекторы в лаборатории?",
        "kb_lab_text": """**Протокол экспериментальной валидации:**

1. **Синтез белка:** рекомбинантная экспрессия в *E. coli* (вектор pET, His-tag)
2. **Очистка:** аффинная хроматография на Ni-NTA колонке
3. **Активность:** измерение термического гистерезиса (наноосмометр Клифтона)
4. **Токсичность:** MTT-тест на клеточных культурах (HEK293, HepG2)
5. **Криоконсервация:** заморозка клеток/тканей с белком (-196°C)
6. **Оценка:** сравнение выживаемости с DMSO-контролем (проточная цитометрия)

**Ключевые показатели успеха:**
- Выживаемость клеток > 80% после разморозки
- Отсутствие кристаллов льда при микроскопии
- Сохранение метаболической активности (АТФ-тест)""",
        "kb_why_title": "🌍 Почему это важно?",
        "kb_why_text": """**Проблема донорства органов:**
- В мире более 150 000 человек находятся в листах ожидания трансплантации
- Каждые 17 минут умирает один пациент, не дождавшийся органа
- До 60% донорских органов не используются из-за логистических ограничений

**Текущие ограничения хранения:**
- Сердце: 4–6 часов
- Почки: 12–24 часа
- Печень: 8–12 часов
- Лёгкие: 4–6 часов

**Потенциал криоконсервации:**
- Хранение органов в течение нескольких лет при -196°C
- Создание «банков органов» по аналогии с банками крови
- Транспортировка органов между континентами
- Спасение тысяч жизней ежегодно""",
        "links_header": "🔗 Полезные ссылки",
        "links_text": "- [UniProt](https://www.uniprot.org/) — база данных белковых последовательностей\n- [PDB](https://www.rcsb.org/) — банк трёхмерных структур белков\n- [ToxinPred2](https://webs.iiitd.edu.in/raghava/toxinpred2/) — предсказание токсичности пептидов\n- [CellPPD](http://crdd.osdd.net/raghava/cellppd/) — предсказание проникающих в клетки пептидов\n- [ESMFold](https://esmatlas.com/) — предсказание 3D-структур белков",
        "footer": "*CryoProtect-Screener v1.0 | Создано для проекта по биотехнологии | 2026*"
    },
    "en": {
        "title": "🧬 CryoProtect-Screener",
        "subtitle": "*In silico* screening of proteins for biomaterial cryopreservation",
        "sidebar_header": "📋 About",
        "sidebar_text": """**Problem:** Every 17 minutes, a patient awaiting transplantation dies. Donor organs can be stored for only 4–24 hours. Cryopreservation at -196°C could extend storage to years, creating "organ banks".

**Solution:** Antifreeze proteins (AFPs), heat shock proteins (HSPs), and LEA proteins — natural cryoprotectants with low toxicity.

**How to use:**
1. Enter a UniProt ID
2. Click "Analyze"
3. Get CryoScore and application recommendations

**Author:** Veronika Bondarenko""",
        "sidebar_footer": "Support science! Organ donation saves lives.",
        "select_preset": "Or choose from examples:",
        "input_label": "UniProt ID:",
        "input_placeholder": "Example: P19479",
        "criteria_header": "🎯 Cryoprotectant Criteria",
        "criteria_text": """**Universal criteria:**
1. ❄️ Ice/water binding
2. 💧 Hydrophilic (GRAVY < 0)
3. 📏 Small size (< 30 kDa)
4. 🧬 Stable (instability index < 40)
5. ✅ Safe (predicted non-toxic)

**Specific criteria:**
- HSP: thermal stability > 60°C
- LEA: high hydrophilicity, desiccation protection
- Peptides: cell membrane penetration""",
        "analyze_button": "🔬 Analyze Protein",
        "warning_id": "⚠️ Please enter a UniProt ID!",
        "loading_sequence": "Loading sequence from UniProt...",
        "error_load": "❌ Failed to load",
        "check_id": "Check the ID.",
        "success_load": "✅ Protein loaded:",
        "type_hsp": "🔬 Protein type: Heat Shock Protein (HSP) — molecular chaperone, stress protection",
        "type_lea": "🔬 Protein type: LEA protein — desiccation protection (anhydrobiosis)",
        "type_afp": "🔬 Protein type: Antifreeze Protein (AFP) — ice binding, thermal hysteresis",
        "type_other": "🔬 Protein type: Synthetic/Other — requires experimental validation",
        "analyzing": "Analyzing physicochemical properties...",
        "error_short": "❌ Sequence too short for analysis.",
        "results_header": "📊 Analysis Results",
        "cryoscore_label": "CryoScore",
        "safety_label": "Safety",
        "penetration_label": "Cell Penetration",
        "excellent_delta": "Excellent candidate!",
        "moderate_delta": "Moderate",
        "low_delta": "Low priority",
        "passport_header": "📋 Cryoprotectant Candidate Passport",
        "passport_param_col": "Parameter",
        "passport_value_col": "Value",
        "passport_score_col": "Assessment",
        "passport_rows": ["UniProt ID", "Length (aa)", "Mol. weight (kDa)", "GRAVY (hydrophobicity)", "T+S+A% (ice-binding)", "Isoelectric point", "Instability index", "Toxicity", "Allergenicity", "In vitro stability", "Cell penetration"],
        "score_small": "✅ Small", "score_large": "⚠️ Large",
        "score_light": "✅ Light", "score_heavy": "⚠️ Heavy",
        "score_hydrophilic": "✅ Hydrophilic", "score_hydrophobic": "⚠️ Hydrophobic",
        "score_high": "✅ High", "score_low": "⚠️ Low",
        "score_stable": "✅ Stable", "score_unstable": "⚠️ Unstable",
        "score_yes": "✅", "score_warn": "⚠️",
        "stable_text": "Stable", "unstable_text": "Unstable",
        "toxic_large": "Non-toxic (large protein)",
        "toxic_short": "Possible toxicity (short peptide)",
        "toxic_likely": "Likely non-toxic",
        "allergen_possible": "Possible allergenicity",
        "allergen_likely": "Likely non-allergenic",
        "safe_yes": "✅ SAFE",
        "safe_check": "⚠️ Requires testing",
        "pen_strong": "✅ CPP+ (strong)",
        "pen_moderate": "🟡 Moderate",
        "pen_weak": "🟡 Weak (small size)",
        "pen_none": "❌ Non-penetrating",
        "cryoscore_chart": "📈 CryoScore Components",
        "comp_ice": "Ice-binding\n(35%)",
        "comp_hydro": "Hydrophilicity\n(25%)",
        "comp_size": "Size\n(20%)",
        "comp_hbond": "H-bonds\n(10%)",
        "comp_stab": "Stability\n(10%)",
        "chart_x": "Component", "chart_y": "Score",
        "chart_cryo_title": "Contribution of Components to CryoScore",
        "kd_header": "💧 Hydrophobicity Profile (Kyte-Doolittle)",
        "kd_x": "Position (window 5 aa)", "kd_y": "Hydrophobicity",
        "kd_title": "Sliding Hydrophobicity (window 5 aa)",
        "kd_boundary": "Hydrophobicity threshold",
        "aa_header": "🧬 Amino Acid Composition",
        "aa_thr": "Threonine (T)", "aa_ser": "Serine (S)", "aa_ala": "Alanine (A)",
        "aa_arg": "Arginine (R)", "aa_lys": "Lysine (K)", "aa_other": "Other",
        "aa_title": "Protein Composition",
        "func_hydrophobic": "Hydrophobic", "func_polar": "Polar (T,S)",
        "func_charged": "Charged (R,K)", "func_other": "Other",
        "func_title": "Functional Groups",
        "aa_col": "Amino acid", "pct_col": "Percentage",
        "group_col": "Group",
        "rec_header": "💡 Recommendation",
        "rec_excellent_title": "🏆 EXCELLENT CANDIDATE!",
        "rec_moderate_title": "🟡 MODERATE CANDIDATE",
        "rec_low_title": "🔴 LOW PRIORITY",
        "rec_moderate_text": "Can be tested, but not a top priority. Consider comparing with other proteins from your library.",
        "rec_low_text": "Consider searching for other proteins. Look for proteins with high T+S+A% and negative GRAVY.",
        "rec_hsp": """**Recommendations for sperm cryopreservation:**
1. Addition to cryoprotective medium (0.1–1 mg/mL)
2. Incubation with spermatozoa (30 min, 37°C)
3. Freezing in liquid nitrogen vapor
4. Post-thaw motility assessment (CASA)""",
        "rec_lea": """**Recommendations for dry preservation:**
1. Electroporation of protein into target cells
2. Desiccation at 60% relative humidity (24 h)
3. Storage at room temperature
4. Rehydration and viability assessment (MTT assay)""",
        "rec_afp": """**Recommendations for organ cryopreservation:**
1. Recombinant protein synthesis in *E. coli*
2. Thermal hysteresis measurement (nanolitre osmometer)
3. MTT assay on HEK293 culture
4. Hepatocyte cryopreservation (0.1–1 mg/mL)
5. Viability comparison with DMSO control""",
        "download_header": "📥 Download Results",
        "download_button": "📄 Download Passport (CSV)",
        "kb_header": "📚 Knowledge Base",
        "kb_afp_title": "ℹ️ What are Antifreeze Proteins (AFPs)?",
        "kb_afp_text": """**Antifreeze proteins (AFPs)** are natural cryoprotectants found in organisms living at subzero temperatures (fish, insects, plants, bacteria).

**Mechanism of action:**
- AFPs adsorb onto the surface of microscopic ice crystals
- They inhibit ice growth (thermal hysteresis effect)
- They lower the freezing point without changing the melting point

**AFP classification:**
- **Type I** (fish): alpha-helical, alanine-rich
- **Type II** (fish): C-type lectins, beta-sandwich
- **Type III** (fish): globular, beta-sandwich
- **AFGP** (fish): glycoproteins with Thr-Ala-Ala repeats
- **Insect:** hyperactive, right-handed beta-helix
- **Plant:** chitinases, thaumatin-like proteins, LTPs
- **Bacterial:** RTX repeats, Ca²⁺-dependent adhesins

**Evolution:** AFPs have arisen convergently at least 6–7 times across different taxa.""",
        "kb_cryo_title": "🧮 How is CryoScore calculated?",
        "kb_cryo_text": """**Formula:** CryoScore = 0.35×TSA + 0.25×(-GRAVY) + 0.20×(1000/MW) + 0.10×HBonds + 0.10×Stability

**Components:**
- **T+S+A% (35%):** ice-binding potential — Thr, Ser, Ala form an ice-complementary surface
- **-GRAVY (25%):** hydrophilicity — more hydrophilic proteins dissolve better in cryoprotective solutions
- **1/MW (20%):** size — low molecular weight proteins penetrate tissues more efficiently
- **H-bonds (10%):** hydrogen bond density correlates with structural stability
- **Stability (10%):** inverted instability index (< 40 = stable in vitro)

**Interpretation:**
- > 70: excellent candidate, recommended for in vitro testing
- 50–70: moderate candidate
- < 50: low priority""",
        "kb_lab_title": "🔬 How are cryoprotectants tested in the lab?",
        "kb_lab_text": """**Experimental validation protocol:**

1. **Protein synthesis:** recombinant expression in *E. coli* (pET vector, His-tag)
2. **Purification:** affinity chromatography on Ni-NTA column
3. **Activity:** thermal hysteresis measurement (Clifton nanolitre osmometer)
4. **Toxicity:** MTT assay on cell cultures (HEK293, HepG2)
5. **Cryopreservation:** freezing of cells/tissues with protein (-196°C)
6. **Evaluation:** viability comparison with DMSO control (flow cytometry)

**Key success indicators:**
- Cell viability > 80% after thawing
- Absence of ice crystals on microscopy
- Preservation of metabolic activity (ATP assay)""",
        "kb_why_title": "🌍 Why does this matter?",
        "kb_why_text": """**The organ donation problem:**
- Over 150,000 people worldwide are on transplant waiting lists
- One patient dies every 17 minutes while waiting for an organ
- Up to 60% of donor organs go unused due to logistical constraints

**Current storage limitations:**
- Heart: 4–6 hours
- Kidneys: 12–24 hours
- Liver: 8–12 hours
- Lungs: 4–6 hours

**Potential of cryopreservation:**
- Storage of organs for years at -196°C
- Creation of "organ banks" similar to blood banks
- Intercontinental organ transportation
- Saving thousands of lives annually""",
        "links_header": "🔗 Useful Links",
        "links_text": "- [UniProt](https://www.uniprot.org/) — protein sequence database\n- [PDB](https://www.rcsb.org/) — 3D protein structure bank\n- [ToxinPred2](https://webs.iiitd.edu.in/raghava/toxinpred2/) — peptide toxicity prediction\n- [CellPPD](http://crdd.osdd.net/raghava/cellppd/) — cell-penetrating peptide prediction\n- [ESMFold](https://esmatlas.com/) — 3D protein structure prediction",
        "footer": "*CryoProtect-Screener v1.0 | Created for a biotechnology project | 2026*"
    }
}

# ============================================================
# LANGUAGE SELECTOR (сначала!)
# ============================================================
if "lang" not in st.session_state:
    st.session_state.lang = "ru"

# Используем пустой контейнер в сайдбаре для переключателя
with st.sidebar:
    lang_choice = st.radio(
        "🌐 Language / Язык",
        ["🇷🇺 Русский", "🇬🇧 English"],
        index=0 if st.session_state.lang == "ru" else 1,
        key="lang_selector"
    )
    st.session_state.lang = "ru" if "Русский" in lang_choice else "en"

T = TEXTS[st.session_state.lang]

# ============================================================
# FUNCTIONS
# ============================================================

@st.cache_data(ttl=3600)
def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200 and response.text.strip():
            record = SeqIO.read(StringIO(response.text), "fasta")
            return str(record.seq), record.description
    except:
        pass
    return None, None


def analyze_properties(sequence):
    clean_seq = sequence.replace("*", "").upper()
    clean_seq = "".join([aa for aa in clean_seq if aa in "ACDEFGHIKLMNPQRSTVWY"])
    
    if len(clean_seq) < 5:
        return None
    
    analysis = ProteinAnalysis(clean_seq)
    
    try:
        aa_percents = analysis.get_amino_acids_percent()
    except AttributeError:
        try:
            aa_percents = analysis.get_amino_acid_percentages()
        except AttributeError:
            counts = Counter(clean_seq)
            total = len(clean_seq)
            aa_percents = {aa: count/total for aa, count in counts.items()}
    
    kd_scale = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
    }
    
    window = 5
    kd_profile = []
    for i in range(len(clean_seq) - window + 1):
        segment = clean_seq[i:i+window]
        score = sum(kd_scale.get(aa, 0) for aa in segment) / window
        kd_profile.append(score)
    
    return {
        "sequence": clean_seq, "length": len(clean_seq),
        "weight": analysis.molecular_weight(), "gravy": analysis.gravy(),
        "instability": analysis.instability_index(),
        "isoelectric_point": analysis.isoelectric_point(),
        "aromaticity": analysis.aromaticity(),
        "thr_pct": aa_percents.get('T', 0) * 100,
        "ser_pct": aa_percents.get('S', 0) * 100,
        "ala_pct": aa_percents.get('A', 0) * 100,
        "tsa_pct": (aa_percents.get('T', 0) + aa_percents.get('S', 0) + aa_percents.get('A', 0)) * 100,
        "arg_pct": aa_percents.get('R', 0) * 100,
        "lys_pct": aa_percents.get('K', 0) * 100,
        "kd_profile": kd_profile,
        "kd_min": min(kd_profile) if kd_profile else 0,
        "kd_max": max(kd_profile) if kd_profile else 0,
    }


def calculate_cryoscore(props):
    w1, w2, w3, w4, w5 = 0.35, 0.25, 0.20, 0.10, 0.10
    
    tsa_score = min(props['tsa_pct'] / 100, 1.0)
    gravy_score = max(-props['gravy'] / 2, 0)
    mw_score = min(1000 / max(props['weight'], 1), 1.0)
    stability_score = max(1.0 - (props['instability'] / 100), 0)
    hbond_score = min(props['tsa_pct'] / 200 + 0.3, 1.0)
    
    cryoscore = (w1 * tsa_score + w2 * gravy_score + w3 * mw_score + 
                 w4 * hbond_score + w5 * stability_score) * 100
    
    return {
        'cryoscore': cryoscore,
        'tsa_score': tsa_score * 100,
        'gravy_score': gravy_score * 100,
        'mw_score': mw_score * 100,
        'hbond_score': hbond_score * 100,
        'stability_score': stability_score * 100
    }


def predict_toxicity(props):
    r = {}
    r['stable_in_vitro'] = ('✅ ' + T['stable_text']) if props['instability'] < 40 else ('⚠️ ' + T['unstable_text'])
    
    if props['length'] > 100:
        r['toxic'] = '✅ ' + T['toxic_large']
    elif props['length'] < 20:
        r['toxic'] = '⚠️ ' + T['toxic_short']
    else:
        r['toxic'] = '✅ ' + T['toxic_likely']
    
    r['allergen'] = ('⚠️ ' + T['allergen_possible']) if props['length'] < 50 else ('✅ ' + T['allergen_likely'])
    
    if '✅' in r['stable_in_vitro'] and '✅' in r['toxic'] and '✅' in r['allergen']:
        r['overall'] = T['safe_yes']
    else:
        r['overall'] = T['safe_check']
    
    return r


def predict_penetration(props):
    pc = (props['arg_pct'] + props['lys_pct']) / 100
    if pc > 0.15 and props['length'] < 200:
        return T['pen_strong']
    elif pc > 0.10:
        return T['pen_moderate']
    elif props['length'] < 100:
        return T['pen_weak']
    else:
        return T['pen_none']


def get_protein_type(uid):
    if uid in ["P0DMV8", "P0A6Y8", "P02829", "P04792"]:
        return T['type_hsp']
    elif uid in ["A0A0K0XZP8", "Q9NFL3"]:
        return T['type_lea']
    elif uid in ["Q9S8C5", "P82972", "A8JHB7"]:
        return T['type_afp']
    else:
        return T['type_other']


def get_recommendation(uid):
    if uid in ["P0DMV8", "P0A6Y8", "P02829", "P04792"]:
        return T['rec_hsp']
    elif uid in ["A0A0K0XZP8", "Q9NFL3"]:
        return T['rec_lea']
    else:
        return T['rec_afp']


# ============================================================
# PRESETS
# ============================================================
PRESETS = {
    "Выберите пример... / Select example...": "",
    "❄️ AFP ржи / Rye AFP (Q9S8C5)": "Q9S8C5",
    "❄️ TmAFP (P82972)": "P82972",
    "❄️ IBP водоросли / Algal IBP (A8JHB7)": "A8JHB7",
    "🔥 HSP70 человека / Human HSP70 (P0DMV8)": "P0DMV8",
    "🔥 HSP70 E. coli (P0A6Y8)": "P0A6Y8",
    "🔥 HSP90 дрожжей / Yeast HSP90 (P02829)": "P02829",
    "🔥 HSP27 человека / Human HSP27 (P04792)": "P04792",
    "💧 LEA тихоходки / Tardigrade LEA (A0A0K0XZP8)": "A0A0K0XZP8",
    "💧 LEA артемии / Artemia LEA (Q9NFL3)": "Q9NFL3",
}

# ============================================================
# HEADER
# ============================================================
st.title(T["title"])
st.markdown(T["subtitle"])

# ============================================================
# SIDEBAR CONTENT (используем контейнеры для динамического обновления)
# ============================================================
with st.sidebar:
    # Разделитель после переключателя языка
    st.markdown("---")
    
    # Заголовок и текст сайдбара
    sidebar_header = st.empty()
    sidebar_text = st.empty()
    st.markdown("---")
    sidebar_footer = st.empty()
    
    # Заполняем актуальными данными
    sidebar_header.header(T["sidebar_header"])
    sidebar_text.markdown(T["sidebar_text"])
    sidebar_footer.markdown(T["sidebar_footer"])

# ============================================================
# MAIN INTERFACE
# ============================================================
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown(f"### 🔍 {T['input_label']}")
    preset = st.selectbox(T["select_preset"], list(PRESETS.keys()))
    uniprot_id = st.text_input(T["input_label"], value=PRESETS.get(preset, ""), placeholder=T["input_placeholder"])

with col2:
    st.markdown(T["criteria_header"])
    st.info(T["criteria_text"])

# ============================================================
# ANALYZE BUTTON & RESULTS
# ============================================================
if st.button(T["analyze_button"], type="primary", use_container_width=True):
    if not uniprot_id:
        st.warning(T["warning_id"])
    else:
        with st.spinner(T["loading_sequence"]):
            sequence, description = fetch_uniprot_sequence(uniprot_id.strip())
        
        if sequence is None:
            st.error(f"{T['error_load']} {uniprot_id}. {T['check_id']}")
        else:
            st.success(f"{T['success_load']} {description[:100]}...")
            st.info(get_protein_type(uniprot_id))
            
            with st.spinner(T["analyzing"]):
                props = analyze_properties(sequence)
            
            if props is None:
                st.error(T["error_short"])
            else:
                scores = calculate_cryoscore(props)
                toxicity = predict_toxicity(props)
                penetration = predict_penetration(props)
                
                # --- RESULTS ---
                st.markdown("---")
                st.markdown(T["results_header"])
                
                col_s, col_f, col_p = st.columns(3)
                
                score = scores['cryoscore']
                with col_s:
                    delta_text = T["excellent_delta"] if score > 70 else (T["moderate_delta"] if score > 50 else T["low_delta"])
                    st.metric(T["cryoscore_label"], f"{score:.1f}/100", delta=delta_text)
                with col_f:
                    st.metric(T["safety_label"], toxicity['overall'])
                with col_p:
                    st.metric(T["penetration_label"], penetration)
                
                # --- PASSPORT ---
                st.markdown(T["passport_header"])
                
                passport_values = [
                    uniprot_id, str(props['length']), f"{props['weight']/1000:.1f}", f"{props['gravy']:.3f}",
                    f"{props['tsa_pct']:.1f}%", f"{props['isoelectric_point']:.2f}", f"{props['instability']:.1f}",
                    toxicity['toxic'], toxicity['allergen'], toxicity['stable_in_vitro'], penetration
                ]
                
                passport_scores = [
                    "—",
                    T["score_small"] if props['length'] < 200 else T["score_large"],
                    T["score_light"] if props['weight'] < 30000 else T["score_heavy"],
                    T["score_hydrophilic"] if props['gravy'] < 0 else T["score_hydrophobic"],
                    T["score_high"] if props['tsa_pct'] > 20 else T["score_low"],
                    "—",
                    T["score_stable"] if props['instability'] < 40 else T["score_unstable"],
                    T["score_yes"] if "✅" in toxicity['toxic'] else T["score_warn"],
                    T["score_yes"] if "✅" in toxicity['allergen'] else T["score_warn"],
                    T["score_yes"] if "✅" in toxicity['stable_in_vitro'] else T["score_warn"],
                    T["score_yes"] if "✅" in penetration else T["score_warn"],
                ]
                
                df = pd.DataFrame({
                    T["passport_param_col"]: T["passport_rows"],
                    T["passport_value_col"]: passport_values,
                    T["passport_score_col"]: passport_scores
                })
                st.dataframe(df, use_container_width=True, hide_index=True)
                
                # --- CRYOSCORE CHART ---
                st.markdown(T["cryoscore_chart"])
                
                comp = {
                    T["comp_ice"]: scores['tsa_score'],
                    T["comp_hydro"]: scores['gravy_score'],
                    T["comp_size"]: scores['mw_score'],
                    T["comp_hbond"]: scores['hbond_score'],
                    T["comp_stab"]: scores['stability_score']
                }
                
                fig = px.bar(x=list(comp.keys()), y=list(comp.values()),
                           labels={'x': T["chart_x"], 'y': T["chart_y"]},
                           title=T["chart_cryo_title"],
                           color=list(comp.values()), color_continuous_scale='RdYlGn', range_color=[0, 100])
                fig.update_layout(showlegend=False, height=400)
                st.plotly_chart(fig, use_container_width=True)
                
                # --- KD PROFILE ---
                st.markdown(T["kd_header"])
                if props['kd_profile']:
                    fig_kd = px.line(y=props['kd_profile'],
                                   labels={'index': T["kd_x"], 'value': T["kd_y"]},
                                   title=T["kd_title"])
                    fig_kd.add_hline(y=0, line_dash="dash", line_color="red", annotation_text=T["kd_boundary"])
                    fig_kd.add_hline(y=props['gravy'], line_dash="dot", line_color="blue",
                                   annotation_text=f"GRAVY={props['gravy']:.3f}")
                    fig_kd.update_layout(height=400)
                    st.plotly_chart(fig_kd, use_container_width=True)
                
                # --- AA COMPOSITION ---
                st.markdown(T["aa_header"])
                col_p1, col_p2 = st.columns(2)
                
                other_pct = 100 - props['thr_pct'] - props['ser_pct'] - props['ala_pct'] - props['arg_pct'] - props['lys_pct']
                
                with col_p1:
                    aa_df = pd.DataFrame({
                        T["aa_col"]: [T["aa_thr"], T["aa_ser"], T["aa_ala"], T["aa_arg"], T["aa_lys"], T["aa_other"]],
                        T["pct_col"]: [props['thr_pct'], props['ser_pct'], props['ala_pct'], props['arg_pct'], props['lys_pct'], other_pct]
                    })
                    fig1 = px.pie(aa_df, values=T["pct_col"], names=T["aa_col"], title=T["aa_title"],
                                 color_discrete_sequence=['#ff6b6b', '#ffd93d', '#6bcb77', '#4d96ff', '#9b59b6', '#bdc3c7'])
                    fig1.update_traces(textposition='inside', textinfo='percent+label')
                    fig1.update_layout(height=400)
                    st.plotly_chart(fig1, use_container_width=True)
                
                with col_p2:
                    func_df = pd.DataFrame({
                        T["group_col"]: [T["func_hydrophobic"], T["func_polar"], T["func_charged"], T["func_other"]],
                        T["pct_col"]: [props['ala_pct'], props['thr_pct'] + props['ser_pct'],
                                      props['arg_pct'] + props['lys_pct'], other_pct]
                    })
                    fig2 = px.pie(func_df, values=T["pct_col"], names=T["group_col"], title=T["func_title"],
                                 color_discrete_sequence=['#e74c3c', '#2ecc71', '#3498db', '#95a5a6'])
                    fig2.update_traces(textposition='inside', textinfo='percent+label')
                    fig2.update_layout(height=400)
                    st.plotly_chart(fig2, use_container_width=True)
                
                # --- RECOMMENDATION ---
                st.markdown("---")
                st.markdown(T["rec_header"])
                
                if score > 70 and '✅' in toxicity['overall']:
                    st.success(T["rec_excellent_title"])
                    st.markdown(get_recommendation(uniprot_id))
                elif score > 50:
                    st.warning(T["rec_moderate_title"])
                    st.markdown(T["rec_moderate_text"])
                else:
                    st.error(T["rec_low_title"])
                    st.markdown(T["rec_low_text"])
                
                # --- DOWNLOAD ---
                st.markdown(T["download_header"])
                csv_data = df.to_csv(index=False).encode('utf-8')
                st.download_button(T["download_button"], csv_data, f"{uniprot_id}_passport.csv", "text/csv")

# ============================================================
# KNOWLEDGE BASE
# ============================================================
st.markdown("---")
st.markdown(T["kb_header"])

with st.expander(T["kb_afp_title"]):
    st.markdown(T["kb_afp_text"])

with st.expander(T["kb_cryo_title"]):
    st.markdown(T["kb_cryo_text"])

with st.expander(T["kb_lab_title"]):
    st.markdown(T["kb_lab_text"])

with st.expander(T["kb_why_title"]):
    st.markdown(T["kb_why_text"])

# ============================================================
# LINKS & FOOTER
# ============================================================
st.markdown("---")
st.markdown(T["links_header"])
st.markdown(T["links_text"])

st.markdown("---")
st.markdown(T["footer"])