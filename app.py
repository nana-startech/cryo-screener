"""
Cold-Adapted Protein Explorer
Educational tool for exploring physicochemical properties of proteins
associated with cold tolerance, with comparison to known AFPs.

DISCLAIMER:
This tool does NOT predict cryoprotective activity. It calculates objective
physicochemical properties and compares them to a reference set of known
antifreeze proteins from public databases.

Author: Veronika Bondarenko
Version: 1.1
"""

import os
import streamlit as st
import requests
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import Counter
import plotly.express as px
import plotly.graph_objects as go
from io import StringIO

st.set_page_config(
    page_title="Cold-Adapted Protein Explorer",
    page_icon="🧊",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================
# REFERENCE DATA — только проверенные значения из ProtParam
# GRAVY рассчитан по реальным последовательностям из UniProt
# TH и IRI — из литературы, где доступно
# ============================================================

KNOWN_AFPS = [
    {"id": "P19479", "name": "Type III AFP", "organism": "Macrozoarces americanus", "length": 65, "mw_kda": 7.0, "gravy": 0.17, "tha_degC": 0.8, "irI": "Moderate"},
    {"id": "P0A3E0", "name": "Type I AFP", "organism": "Pseudopleuronectes americanus", "length": 37, "mw_kda": 4.5, "gravy": 0.50, "tha_degC": 0.7, "irI": "Low"},
    {"id": "P82972", "name": "Hyperactive AFP", "organism": "Tenebrio molitor", "length": 70, "mw_kda": 7.5, "gravy": -0.32, "tha_degC": 5.5, "irI": "High"},
    {"id": "A8JHB7", "name": "Ice-binding protein", "organism": "Chlamydomonas sp.", "length": 250, "mw_kda": 27.5, "gravy": -0.41, "tha_degC": 2.1, "irI": "Moderate"},
    {"id": "P0DMV8", "name": "HSP70 (control)", "organism": "Homo sapiens", "length": 641, "mw_kda": 70.0, "gravy": -0.45, "tha_degC": None, "irI": None},
    {"id": "P0A6Y8", "name": "HSP70 (control)", "organism": "Escherichia coli", "length": 638, "mw_kda": 69.0, "gravy": -0.43, "tha_degC": None, "irI": None},
    {"id": "A0A0K0XZP8", "name": "LEA protein", "organism": "Tardigrada", "length": 200, "mw_kda": 22.0, "gravy": -1.20, "tha_degC": None, "irI": None},
]

# ============================================================
# TRANSLATIONS
# ============================================================

TEXTS = {
    "ru": {
        "title": "🧬 Cold-Adapted Protein Explorer",
        "subtitle": "Анализ физико-химических свойств белков и сравнение с известными антифризными белками",
        "main_warning": """### ⚠️ Важно: это образовательный инструмент, а не предиктор AFP

Этот сайт **НЕ определяет**, является ли белок антифризным.

Аминокислотный состав и физико-химические параметры **не являются критерием AFP**.

Активность AFP определяется **только экспериментально**:
- **Термический гистерезис (TH)** — разница между замерзанием и плавлением
- **Ингибирование рекристаллизации льда (IRI)** — способность предотвращать рост крупных кристаллов""",
        "sidebar_header": "📋 О проекте",
        "sidebar_text": """**Что делает:**
- Загружает последовательность из UniProt
- Вычисляет физико-химические свойства
- Показывает их рядом с известными AFP
- Предсказывает структуру (с оговорками)

**Чего не делает:**
- Не предсказывает активность
- Не определяет AFP
- Не заменяет лабораторные тесты

**Автор:** Бондаренко Вероника
**Версия:** 3.1""",
        "sidebar_footer": "Образовательный прототип. Не для клинического использования.",
        "select_preset": "Выберите пример:",
        "input_label": "UniProt ID:",
        "input_placeholder": "Например: P19479",
        "analyze_button": "🔬 Анализировать",
        "warning_id": "⚠️ Введите UniProt ID!",
        "loading_sequence": "Загрузка из UniProt...",
        "error_load": "❌ Ошибка загрузки",
        "check_id": "Проверьте ID.",
        "success_load": "✅ Загружен:",
        "analyzing": "Анализ...",
        "error_short": "❌ Последовательность слишком короткая.",
        "nonstandard_aa_warning": "⚠️ Обнаружены нестандартные аминокислоты",
        "nonstandard_aa_note": "Следующие символы были заменены: B→N, Z→Q, X→A, J→L, U→C, O→K. Это стандартная практика в биоинформатике. Результаты могут быть приблизительными.",
        "results_header": "📊 Физико-химические свойства",
        "passport_header": "📋 Свойства белка",
        "passport_param_col": "Параметр",
        "passport_value_col": "Значение",
        "passport_note_col": "Комментарий",
        "passport_rows": [
            "UniProt ID", "Длина (а.о.)", "Мол. масса (кДа)",
            "GRAVY", "pI", "Индекс нестабильности",
            "Алифатический индекс", "Заряд при pH 7.0",
            "Cys%", "Ароматичность"
        ],
        "passport_notes": [
            "—",
            "Объективный параметр",
            "Объективный параметр",
            "Средняя гидрофобность. AFP могут быть как гидрофобными, так и гидрофильными",
            "Объективный параметр",
            "Предсказывает стабильность при 37°C, не при заморозке",
            "Показатель термостабильности",
            "Упрощённый расчёт. Не учитывает окружение белка",
            "Важен для дисульфидных связей",
            "Относительное содержание ароматических аминокислот"
        ],
        "comparison_header": "🔄 Сравнение с известными AFP",
        "comparison_text": "Значения вашего белка показаны вместе с известными AFP. Это сравнение **НЕ означает**, что ваш белок является антифризным. Оно лишь показывает, попадают ли численные значения в диапазон известных AFP.",
        "ref_table_header": "📊 Референсная таблица",
        "ref_table_note": "Примечание: значения TH и IRI приведены не для всех белков из-за неполноты экспериментальных данных в открытой литературе.",
        "your_protein_label": "Ваш белок",
        "known_afp_label": "Известные AFP / контроли",
        "grav_dist_header": "📊 Распределение по GRAVY",
        "grav_dist_text": "Гистограмма значений GRAVY для известных AFP. Ваш белок показан красной линией.",
        "kd_header": "💧 Профиль гидрофобности (Kyte-Doolittle)",
        "kd_x": "Позиция (окно 5 а.о.)",
        "kd_y": "Гидрофобность",
        "kd_title": "Скользящая гидрофобность",
        "kd_boundary": "Граница",
        "aa_header": "🧬 Аминокислотный состав",
        "aa_col": "Аминокислота",
        "pct_col": "Процент",
        "aa_title": "Состав белка",
        "structure_header": "🧠 Предсказанная структура",
        "structure_note": "Структура предсказана через ESMFold API. Файл доступен, если средний pLDDT > 70. **Важно:** pLDDT > 70 не является абсолютной гарантией правильности структуры.",
        "structure_button": "Скачать PDB",
        "structure_error": "Качественное предсказание структуры недоступно для этого белка.",
        "structure_plddt_label": "Средний pLDDT:",
        "download_header": "📥 Скачать",
        "download_button": "📄 Скачать CSV",
        "kb_header": "📚 Научная справка",
        "kb_th_title": "🌡️ Что такое термический гистерезис (TH)?",
        "kb_th_text": """**Термический гистерезис (TH)** — разница между точкой замерзания и точкой плавления раствора.

В присутствии AFP лёд начинает расти только при температуре значительно ниже точки плавления. Эта разница называется TH.

**Пример:** AFP Type III из Ocean pout даёт TH ~0.7–0.8°C при ~10 мг/мл.

**Важно:** TH и IRI — разные активности.""",
        "kb_iri_title": "🧊 Что такое ингибирование рекристаллизации (IRI)?",
        "kb_iri_text": """**Ингибирование рекристаллизации льда (IRI)** — способность белка предотвращать рост крупных кристаллов льда за счёт мелких.

При заморозке-оттаивании мелкие кристаллы сливаются в крупные, которые повреждают клетки. AFP подавляют этот процесс.

**Важно:** белок может иметь сильную IRI-активность при слабом TH, и наоборот.

**Эти активности измеряются только экспериментально.**""",
        "limitations_header": "⚠️ Ограничения",
        "limitations_text": """Этот инструмент — образовательный. Он не предсказывает активность.

Физико-химические параметры **не коррелируют напрямую** с криопротекторной активностью.

Единственный способ определить, является ли белок AFP — провести экспериментальные измерения TH и IRI.""",
        "links_header": "🔗 Полезные ссылки",
        "links_text": "- [UniProt](https://www.uniprot.org/)\n- [ESMFold](https://esmatlas.com/)\n- [ToxinPred2](https://webs.iiitd.edu.in/raghava/toxinpred2/)",
        "footer": "*Cold-Adapted Protein Explorer v3.1 | Образовательный прототип | 2026*"
    },
    "en": {
        "title": "🧬 Cold-Adapted Protein Explorer",
        "subtitle": "Analysis of physicochemical properties and comparison with known antifreeze proteins",
        "main_warning": """### ⚠️ Important: this is an educational tool, not an AFP predictor

This site does **NOT determine** whether a protein is antifreeze.

Amino acid composition and physicochemical parameters **are not criteria for AFP**.

AFP activity is determined **only experimentally**:
- **Thermal hysteresis (TH)** — difference between freezing and melting
- **Ice recrystallization inhibition (IRI)** — ability to prevent large crystal growth""",
        "sidebar_header": "📋 About",
        "sidebar_text": """**What it does:**
- Loads sequence from UniProt
- Calculates physicochemical properties
- Shows them alongside known AFPs
- Predicts structure (with caveats)

**What it does NOT do:**
- Does not predict activity
- Does not determine AFP
- Does not replace lab tests

**Author:** Veronika Bondarenko
**Version:** 3.1""",
        "sidebar_footer": "Educational prototype. Not for clinical use.",
        "select_preset": "Choose an example:",
        "input_label": "UniProt ID:",
        "input_placeholder": "Example: P19479",
        "analyze_button": "🔬 Analyze",
        "warning_id": "⚠️ Enter a UniProt ID!",
        "loading_sequence": "Loading from UniProt...",
        "error_load": "❌ Load error",
        "check_id": "Check ID.",
        "success_load": "✅ Loaded:",
        "analyzing": "Analyzing...",
        "error_short": "❌ Sequence too short.",
        "nonstandard_aa_warning": "⚠️ Non-standard amino acids detected",
        "nonstandard_aa_note": "The following characters were replaced: B→N, Z→Q, X→A, J→L, U→C, O→K. This is standard bioinformatics practice. Results may be approximate.",
        "results_header": "📊 Physicochemical Properties",
        "passport_header": "📋 Protein Properties",
        "passport_param_col": "Parameter",
        "passport_value_col": "Value",
        "passport_note_col": "Note",
        "passport_rows": [
            "UniProt ID", "Length (aa)", "Mol. weight (kDa)",
            "GRAVY", "pI", "Instability index",
            "Aliphatic index", "Charge at pH 7.0",
            "Cys%", "Aromaticity"
        ],
        "passport_notes": [
            "—",
            "Objective",
            "Objective",
            "Average hydrophobicity. AFPs can be either hydrophobic or hydrophilic",
            "Objective",
            "Predicts stability at 37°C, not under freezing",
            "Thermostability indicator",
            "Simplified calculation. Does not account for protein environment",
            "Important for disulfide bonds",
            "Relative aromatic content"
        ],
        "comparison_header": "🔄 Comparison with Known AFPs",
        "comparison_text": "Your protein's values are shown alongside known AFPs. This comparison does **NOT mean** your protein is antifreeze. It only shows whether numeric values fall within the range of known AFPs.",
        "ref_table_header": "📊 Reference Table",
        "ref_table_note": "Note: TH and IRI values are not available for all proteins due to incomplete experimental data in open literature.",
        "your_protein_label": "Your protein",
        "known_afp_label": "Known AFPs / controls",
        "grav_dist_header": "📊 GRAVY Distribution",
        "grav_dist_text": "Histogram of GRAVY values for known AFPs. Your protein is shown as a red line.",
        "kd_header": "💧 Hydrophobicity Profile (Kyte-Doolittle)",
        "kd_x": "Position (window 5 aa)",
        "kd_y": "Hydrophobicity",
        "kd_title": "Sliding Hydrophobicity",
        "kd_boundary": "Threshold",
        "aa_header": "🧬 Amino Acid Composition",
        "aa_col": "Amino acid",
        "pct_col": "Percentage",
        "aa_title": "Protein Composition",
        "structure_header": "🧠 Predicted Structure",
        "structure_note": "Structure predicted via ESMFold API. File is available if average pLDDT > 70. **Important:** pLDDT > 70 is not an absolute guarantee of structure correctness.",
        "structure_button": "Download PDB",
        "structure_error": "A reliable structure prediction is unavailable for this protein.",
        "structure_plddt_label": "Average pLDDT:",
        "download_header": "📥 Download",
        "download_button": "📄 Download CSV",
        "kb_header": "📚 Scientific Notes",
        "kb_th_title": "🌡️ What is Thermal Hysteresis (TH)?",
        "kb_th_text": """**Thermal hysteresis (TH)** is the difference between the freezing point and melting point of a solution.

In the presence of AFPs, ice starts growing only at a temperature significantly below the melting point. This difference is called TH.

**Example:** AFP Type III from Ocean pout gives TH ~0.7–0.8°C at ~10 mg/mL.

**Important:** TH and IRI are different activities.""",
        "kb_iri_title": "🧊 What is Ice Recrystallization Inhibition (IRI)?",
        "kb_iri_text": """**Ice recrystallization inhibition (IRI)** is the ability of a protein to prevent large ice crystals from growing at the expense of small ones.

During freeze-thaw, small crystals merge into large ones, which damage cells. AFPs suppress this process.

**Important:** a protein can have strong IRI activity with weak TH, and vice versa.

**These activities are measured only experimentally.**""",
        "limitations_header": "⚠️ Limitations",
        "limitations_text": """This tool is educational. It does not predict activity.

Physicochemical parameters **do not directly correlate** with cryoprotective activity.

The only way to determine whether a protein is AFP is to perform experimental TH and IRI measurements.""",
        "links_header": "🔗 Useful Links",
        "links_text": "- [UniProt](https://www.uniprot.org/)\n- [ESMFold](https://esmatlas.com/)\n- [ToxinPred2](https://webs.iiitd.edu.in/raghava/toxinpred2/)",
        "footer": "*Cold-Adapted Protein Explorer v3.1 | Educational prototype | 2026*"
    }
}

# ============================================================
# LANGUAGE SELECTOR
# ============================================================
if "lang" not in st.session_state:
    st.session_state.lang = "ru"

with st.sidebar:
    lang_choice = st.radio(
        "🌐 Language / Язык",
        ["🇷🇺 Русский", "en English"],
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
    filename = f"{uniprot_id}_cold.fasta"
    if os.path.exists(filename):
        try:
            record = SeqIO.read(filename, "fasta")
            return str(record.seq), record.description
        except:
            pass

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        response = requests.get(url, timeout=15)
        if response.status_code == 200 and response.text.strip():
            record = SeqIO.read(StringIO(response.text), "fasta")
            return str(record.seq), record.description
    except:
        pass

    return None, None


def compute_aliphatic_index(sequence):
    length = len(sequence)
    if length == 0:
        return 0
    counts = Counter(sequence)
    ai = (counts.get('A', 0) * 2.9 +
          counts.get('V', 0) * 3.9 +
          counts.get('L', 0) * 3.9 +
          counts.get('I', 0) * 3.9) / length
    return round(ai * 100, 2)


def compute_charge_at_ph7(sequence):
    charge = 0
    for aa in sequence:
        if aa == 'K':
            charge += 1 / (1 + 10**(7.0 - 10.5))
        elif aa == 'R':
            charge += 1 / (1 + 10**(7.0 - 12.5))
        elif aa == 'H':
            charge += 1 / (1 + 10**(7.0 - 6.0))
        elif aa == 'D':
            charge -= 1 / (1 + 10**(3.9 - 7.0))
        elif aa == 'E':
            charge -= 1 / (1 + 10**(4.1 - 7.0))
        elif aa == 'C':
            charge -= 1 / (1 + 10**(8.3 - 7.0))
        elif aa == 'Y':
            charge -= 1 / (1 + 10**(10.1 - 7.0))
    return round(charge, 3)


def analyze_properties(sequence):
    raw_seq = sequence.replace("*", "").upper()
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    nonstandard = set(raw_seq) - standard_aa

    clean_seq = raw_seq.replace('B', 'N').replace('Z', 'Q').replace('X', 'A').replace('J', 'L').replace('U', 'C').replace('O', 'K')

    if len(clean_seq) < 5:
        return None, nonstandard

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

    props = {
        "sequence": clean_seq,
        "length": len(clean_seq),
        "weight": analysis.molecular_weight(),
        "gravy": analysis.gravy(),
        "instability": analysis.instability_index(),
        "isoelectric_point": analysis.isoelectric_point(),
        "aromaticity": analysis.aromaticity(),
        "aliphatic_index": compute_aliphatic_index(clean_seq),
        "charge_ph7": compute_charge_at_ph7(clean_seq),
        "cys_pct": aa_percents.get('C', 0) * 100,
        "kd_profile": kd_profile,
    }

    return props, nonstandard


def get_structure_with_plddt(uniprot_id):
    url = f"https://api.esmatlas.com/fetchPredictedStructure/{uniprot_id}"
    try:
        response = requests.get(url, timeout=20)
        if response.status_code == 200 and len(response.text) > 100:
            text = response.text
            plddt_values = []
            for line in text.splitlines():
                if line.startswith("ATOM"):
                    b_factor = float(line[60:66].strip())
                    plddt_values.append(b_factor)
            avg_plddt = sum(plddt_values) / len(plddt_values) if plddt_values else 0
            return text, round(avg_plddt, 2)
    except:
        pass
    return None, None


# ============================================================
# PRESETS
# ============================================================
PRESETS = {
    "Выберите пример... / Select example...": "",
    "P19479 (Type III AFP)": "P19479",
    "P0A3E0 (Type I AFP)": "P0A3E0",
    "P82972 (TmAFP)": "P82972",
    "A8JHB7 (IBP водоросли)": "A8JHB7",
    "P0DMV8 (HSP70 человека)": "P0DMV8",
    "P0A6Y8 (HSP70 E. coli)": "P0A6Y8",
    "A0A0K0XZP8 (LEA тихоходки)": "A0A0K0XZP8",
}

# ============================================================
# UI
# ============================================================
st.title(T["title"])
st.markdown(T["subtitle"])
st.markdown(T["main_warning"])

with st.sidebar:
    st.markdown("---")
    st.header(T["sidebar_header"])
    st.markdown(T["sidebar_text"])
    st.markdown("---")
    st.markdown(T["sidebar_footer"])

col1, col2 = st.columns([2, 1])

with col1:
    st.markdown(f"### 🔍 {T['input_label']}")
    preset = st.selectbox(T["select_preset"], list(PRESETS.keys()))
    uniprot_id = st.text_input(T["input_label"], value=PRESETS.get(preset, ""), placeholder=T["input_placeholder"])

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

            with st.spinner(T["analyzing"]):
                props, nonstandard = analyze_properties(sequence)

            if props is None:
                st.error(T["error_short"])
            else:
                if nonstandard:
                    st.warning(T["nonstandard_aa_warning"])
                    st.info(T["nonstandard_aa_note"])

                st.markdown("---")
                st.markdown(T["results_header"])

                values = [
                    uniprot_id,
                    str(props['length']),
                    f"{props['weight']/1000:.1f}",
                    f"{props['gravy']:.3f}",
                    f"{props['isoelectric_point']:.2f}",
                    f"{props['instability']:.1f}",
                    f"{props['aliphatic_index']:.2f}",
                    f"{props['charge_ph7']:.2f}",
                    f"{props['cys_pct']:.1f}%",
                    f"{props['aromaticity']:.3f}"
                ]

                df = pd.DataFrame({
                    T["passport_param_col"]: T["passport_rows"],
                    T["passport_value_col"]: values,
                    T["passport_note_col"]: T["passport_notes"]
                })
                st.dataframe(df, use_container_width=True, hide_index=True)

                st.markdown(T["comparison_header"])
                st.markdown(T["comparison_text"])

                ref_df = pd.DataFrame(KNOWN_AFPS)
                st.markdown(T["ref_table_header"])
                st.markdown(T["ref_table_note"])
                st.dataframe(ref_df[["id", "name", "organism", "length", "mw_kda", "gravy", "tha_degC", "irI"]], use_container_width=True, hide_index=True)

                # GRAVY Distribution
                st.markdown(T["grav_dist_header"])
                st.markdown(T["grav_dist_text"])
                fig_grav = go.Figure()
                fig_grav.add_trace(go.Histogram(x=ref_df['gravy'], nbinsx=8, marker_color='lightblue', opacity=0.7))
                fig_grav.add_vline(x=props['gravy'], line_width=3, line_color="red", annotation_text=f"{props['gravy']:.3f}")
                fig_grav.update_layout(height=300, xaxis_title="GRAVY", yaxis_title="Count")
                st.plotly_chart(fig_grav, use_container_width=True)

                # KD Profile
                st.markdown(T["kd_header"])
                if props['kd_profile']:
                    fig_kd = px.line(
                        y=props['kd_profile'],
                        labels={'index': T["kd_x"], 'value': T["kd_y"]},
                        title=T["kd_title"]
                    )
                    fig_kd.add_hline(y=0, line_dash="dash", line_color="red", annotation_text=T["kd_boundary"])
                    fig_kd.update_layout(height=400)
                    st.plotly_chart(fig_kd, use_container_width=True)

                # AA Composition
                st.markdown(T["aa_header"])
                counts = Counter(props["sequence"])
                aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
                aa_pcts = {aa: counts.get(aa, 0) / len(props["sequence"]) * 100 for aa in aa_order if counts.get(aa, 0) > 0}
                aa_df = pd.DataFrame({
                    T["aa_col"]: list(aa_pcts.keys()),
                    T["pct_col"]: list(aa_pcts.values())
                })
                fig_aa = px.pie(aa_df, values=T["pct_col"], names=T["aa_col"], title=T["aa_title"])
                fig_aa.update_traces(textposition='inside', textinfo='percent+label')
                st.plotly_chart(fig_aa, use_container_width=True)

                # Structure
                st.markdown(T["structure_header"])
                st.markdown(T["structure_note"])
                structure, avg_plddt = get_structure_with_plddt(uniprot_id)
                if structure and avg_plddt > 70:
                    st.success(f"{T['structure_plddt_label']} {avg_plddt}")
                    st.download_button(T["structure_button"], structure, file_name=f"{uniprot_id}.pdb", mime="chemical/x-pdb")
                else:
                    st.info(T["structure_error"])

                # Download CSV
                st.markdown(T["download_header"])
                csv_data = df.to_csv(index=False).encode('utf-8')
                st.download_button(T["download_button"], csv_data, f"{uniprot_id}_properties.csv", "text/csv")

st.markdown("---")
st.markdown(T["kb_header"])

with st.expander(T["kb_th_title"]):
    st.markdown(T["kb_th_text"])

with st.expander(T["kb_iri_title"]):
    st.markdown(T["kb_iri_text"])

with st.expander(T["limitations_header"]):
    st.markdown(T["limitations_text"])

st.markdown("---")
st.markdown(T["links_header"])
st.markdown(T["links_text"])

st.markdown("---")
st.markdown(T["footer"])
