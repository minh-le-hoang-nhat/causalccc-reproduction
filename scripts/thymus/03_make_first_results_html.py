"""Build results/thymus/first_results.html: fig. S11 reproduction, CellPhoneDB by stage, overview, Table S5.
Needs the outputs of 02_run_cellphonedb.py, 04_cpdb_by_stage.py and 05_figS11_expression.py.
Run in WSL: python 03_make_first_results_html.py"""
import base64, html, io
import pandas as pd

R = '/mnt/c/Users/Dell 5540/Documents/computational_biology_project/causalccc-reproduction'
RES = f'{R}/results/thymus'
STAGES = ['early', 'middle', 'late']
STROMA = {'Fb_1', 'Fb_2', 'Fb_cycling', 'Endo', 'VSMC', 'Lymph', 'cTEC', 'mTEC(I)', 'mTEC(II)', 'mTEC(III)',
          'mTEC(IV)', 'mcTEC', 'TEC(neuro)', 'TEC(myo)', 'Epi_GCM2'}
NON_IMMUNE = STROMA | {'Ery', 'Mgk'}
esc = html.escape
split = lambda g: g.rsplit('_', 1)            # 'mTEC(II)_early' -> ['mTEC(II)', 'early']


def fmt_groups(df, col='frac'):
    return ', '.join(f'{esc(g)} ({v:.0%})' for g, v in zip(df['group'], df[col]))


# ---------------------------------------------------------------- step 1: fig. S11
E = pd.read_csv(f'{RES}/figS11/expression_by_type_stage.tsv', sep='\t')
NT = pd.read_csv(f'{RES}/figS11/notch_receptors_ETP_DN.tsv', sep='\t')
g = lambda gene: E[E.gene == gene].sort_values('mean', ascending=False)

ltb_top = g('LTB').head(10)
ltb_immune = (~ltb_top.cell_type.isin(NON_IMMUNE)).sum()
ltbr = g('LTBR')
ltbr_str = ltbr[ltbr.cell_type.isin(STROMA)]
ltbr_imm = ltbr[~ltbr.cell_type.isin(NON_IMMUNE)]
c1 = (f'LTB : {ltb_immune}/10 des groupes les plus exprimants sont immuns ({fmt_groups(ltb_top.head(4))}…). '
      f'LTBR détecté dans ≥20 % des cellules pour {(ltbr_str.frac >= .2).sum()}/{len(ltbr_str)} groupes stromaux '
      f'contre {(ltbr_imm.frac >= .2).sum()}/{len(ltbr_imm)} groupes immuns (surtout DC, macrophages, monocytes).')
v1 = ltb_immune == 10 and (ltbr_str.frac >= .2).mean() > .5

rankl = g('TNFSF11'); rankl20 = rankl[rankl.frac >= .2]
rank = g('TNFRSF11A'); rank20 = rank[rank.frac >= .2]
c2 = (f'TNFSF11 (RANKL) dans ≥20 % des cellules : {fmt_groups(rankl20)} uniquement. '
      f'TNFRSF11A (RANK) dans ≥20 % : {fmt_groups(rank20)}.')
v2 = set(rankl20.cell_type) == {'ILC3'} and {'mTEC(II)', 'Lymph'} <= set(rank20.cell_type)

fgf7 = g('FGF7'); fgf7_50 = fgf7[fgf7.frac >= .5]
fgfr2 = E[(E.gene == 'FGFR2') & E.cell_type.isin(['cTEC', 'mTEC(I)', 'mTEC(II)', 'mTEC(III)', 'mcTEC'])]
fgfr2_tab = fgfr2.pivot(index='cell_type', columns='stage', values='mean').reindex(columns=STAGES)
fgfr2_txt = '; '.join(f'{esc(t)} ' + ' → '.join('—' if pd.isna(fgfr2_tab.loc[t, s]) else f'{fgfr2_tab.loc[t, s]:.2f}'
                                                   for s in STAGES) for t in fgfr2_tab.index)
c3 = (f'FGF7 dans ≥50 % des cellules : {", ".join(sorted(set(fgf7_50.cell_type)))}. '
      f'FGFR2, expression moyenne early → middle → late : {fgfr2_txt}.')
m1 = fgfr2_tab.loc['mTEC(I)']
v3 = set(fgf7_50.cell_type) <= {'Fb_1', 'Fb_2', 'Fb_cycling', 'VSMC'} and m1['late'] < m1['middle']

etp = NT[NT.group.str.startswith('ETP_')].pivot(index='group', columns='gene', values='mean')
notch_ok = all(etp.loc[r].idxmax() == 'NOTCH1' for r in etp.index)
etp_txt = '; '.join(f'{esc(r)} : NOTCH1 {etp.loc[r, "NOTCH1"]:.2f} vs NOTCH2–4 ≤ '
                    f'{etp.loc[r, ["NOTCH2", "NOTCH3", "NOTCH4"]].max():.2f}' for r in etp.index)
jag2_top, dll4_top = g('JAG2').head(3), g('DLL4').head(3)
jag1 = g('JAG1'); jag1_35 = jag1[jag1.frac >= .35]
c4 = (f'Dans les ETP, NOTCH1 est le récepteur Notch le plus exprimé à chaque stade ({etp_txt}). '
      f'Top JAG2 : {fmt_groups(jag2_top)}. Top DLL4 : {fmt_groups(dll4_top)}. '
      f'JAG1 dans ≥35 % des cellules : {fmt_groups(jag1_35)}.')
v4 = notch_ok and {'cTEC', 'Endo'} <= set(jag2_top.cell_type) | set(dll4_top.cell_type)

CLAIMS = [
    ('LTB:LTBR', 'LTB vient de cellules immunes variées ; LTBR est reçu par la plupart des cellules stromales.', c1, v1,
     'LTBR est aussi exprimé par les DC et les macrophages.'),
    ('RANKL:RANK', 'Signal limité aux ILC3 → mTEC(II) / endothélium lymphatique.', c2, v2, ''),
    ('FGF7:FGFR2', 'FGF7 vient des fibroblastes vers les TEC ; FGFR2 baisse chez l\'adulte.', c3, v3,
     'FGF7 aussi dans les VSMC. Baisse testable sur les mTEC seulement : les cTEC late ont <50 cellules.'),
    ('Notch', 'NOTCH1 = principal récepteur des ETP ; JAG2 et DLL4 sur cTEC et endothélium ; JAG1 sur les autres TEC.', c4, v4,
     'JAG1 est aussi fort dans les VSMC.'),
]
claim_rows = ''.join(
    f'<tr><td><b>{esc(a)}</b></td><td>{esc(b)}</td><td>{c}</td>'
    f'<td class="{"ok" if v else "ko"}">{"retrouvé" if v else "à discuter"}</td><td>{esc(n)}</td></tr>'
    for a, b, c, v, n in CLAIMS)
png = base64.b64encode(open(f'{RES}/figS11/figS11_dotplot.png', 'rb').read()).decode()

# ---------------------------------------------------------------- step 2: CellPhoneDB by stage
PS = None
for st in STAGES:   # one CellPhoneDB run per stage, merged on the interaction id
    t = pd.read_csv(f'{RES}/cpdb_by_stage/{st}/statistical_analysis_pvalues_{st}.txt', sep='\t')
    t = t[['id_cp_interaction', 'interacting_pair'] + [c for c in t.columns if '|' in c]]
    PS = t if PS is None else PS.merge(t.drop(columns='interacting_pair'), on='id_cp_interaction')
ps_cols = [c for c in PS.columns if '|' in c]
ps_names = PS['interacting_pair'].astype(str)


def pv(pair, sender, receiver, stage):
    col = f'{sender}_{stage}|{receiver}_{stage}'
    if col not in PS.columns:
        return '<td class="na">absent</td>'
    p = float(PS.loc[ps_names == pair, col].iloc[0])
    return f'<td class="{"ok" if p < .05 else "ko"}">{p:.3f}</td>'


def axis_table(pair, couples):
    h = (f'<h3><code>{esc(pair)}</code></h3><div class="scroll"><table><thead><tr><th>sender → receiver</th>'
         + ''.join(f'<th>{s}</th>' for s in STAGES) + '</tr></thead><tbody>')
    for s, r in couples:
        h += f'<tr><td>{esc(s)} → {esc(r)}</td>' + ''.join(pv(pair, s, r, st) for st in STAGES) + '</tr>'
    return h + '</tbody></table></div>'


TECS = ['cTEC', 'mTEC(I)', 'mTEC(II)', 'mTEC(III)', 'mcTEC']
step2 = axis_table('TNFSF11_TNFRSF11A', [('ILC3', 'mTEC(II)'), ('ILC3', 'Lymph')])
step2 += axis_table('FGF7_FGFR2', [(f, t) for f in ('Fb_1', 'Fb_2') for t in TECS])
for lig in ('DLL4', 'JAG2', 'JAG1'):
    step2 += axis_table(f'{lig}_NOTCH1', [(s, 'ETP') for s in ['cTEC', 'Endo'] + TECS[1:]])

# LTB_LTBR: share of significant same-stage couples that go immune -> stroma
row = PS.loc[ps_names == 'LTB_LTBR', ps_cols].iloc[0]
lt = []
for st in STAGES:
    cs = [c for c in ps_cols if all(split(x)[1] == st for x in c.split('|')) and row[c] < .05]
    i2s = [c for c in cs if split(c.split('|')[0])[0] not in NON_IMMUNE and split(c.split('|')[1])[0] in STROMA]
    lt.append(f'<tr><td>{st}</td><td>{len(cs)}</td><td>{len(i2s)} ({len(i2s) / max(len(cs), 1):.0%})</td></tr>')
step2 += ('<h3><code>LTB_LTBR</code></h3><div class="scroll"><table><thead><tr><th>stade</th>'
          '<th>couples significatifs (même stade)</th><th>dont cellule immune → stroma</th></tr></thead><tbody>'
          + ''.join(lt) + '</tbody></table></div>')

# ---------------------------------------------------------------- overview (all stages pooled), oriented L-R
P = pd.read_csv(f'{RES}/cpdb/statistical_analysis_pvalues_thymus.txt', sep='\t')
pc = [c for c in P.columns if '|' in c]
LR = P[P.directionality == 'Ligand-Receptor']
SENDERS = ['cTEC', 'mTEC(I)', 'mTEC(II)', 'mcTEC', 'DC1', 'DC2', 'aDC', 'pDC', 'Fb_1', 'Fb_2', 'Endo']
RECEIVERS = ['DN', 'DP', 'αβT(entry)', 'CD4+T', 'CD8+T', 'CD4+Tmem', 'CD8+Tmem', 'Treg', 'T(agonist)', 'CD8αα', 'γδT']
lig_first = ~LR.receptor_a.astype(str).eq('True')


def n_lr(s, r):   # interactions with the ligand in s and the receptor in r
    return int(((LR[f'{s}|{r}'] < .05) & lig_first).sum() + ((LR[f'{r}|{s}'] < .05) & ~lig_first).sum())


mat = pd.DataFrame({r: {s: n_lr(s, r) for s in SENDERS} for r in RECEIVERS})
vmax = mat.values.max()
heat = '<table class="heat"><thead><tr><th>ligand ↓ / récepteur →</th>' + ''.join(f'<th>{esc(r)}</th>' for r in RECEIVERS) + '</tr></thead><tbody>'
for s in SENDERS:
    heat += f'<tr><th>{esc(s)}</th>' + ''.join(
        f'<td style="background:rgba(31,111,235,{0.08 + 0.85 * mat.loc[s, r] / vmax:.2f})">{mat.loc[s, r]}</td>' for r in RECEIVERS) + '</tr>'
heat += '</tbody></table>'

# ---------------------------------------------------------------- Table S5, compared by partner ids (order-free)
raw = open(f'{R}/data/thymus/aay3224_table-s5.csv', encoding='utf-8', errors='replace', newline='').read()
raw = raw.replace('\r\n', '\n').replace('\r', '\n')
S5 = pd.read_csv(io.StringIO('\n'.join(raw.split('\n')[1:])), index_col=0)
key = lambda a, b: frozenset([str(a).replace('simple:', ''), str(b).replace('simple:', '')])
s5k = {key(a, b) for a, b in zip(S5.partner_a, S5.partner_b)}
ours = {key(a, b) for a, b in zip(P.partner_a, P.partner_b)}
sig_any = (P[pc] < .05).any(axis=1)
oursig = {key(a, b) for a, b, s in zip(P.partner_a, P.partner_b, sig_any) if s}

CSS = """
:root{--bg:#fff;--fg:#1f2328;--mut:#656d76;--line:#d0d7de;--card:#f6f8fa;--acc:#1f6feb;--ok:#1a7f37;--okbg:#dafbe1;--ko:#9a6700;--kobg:#fff8c5}
@media(prefers-color-scheme:dark){:root{--bg:#0d1117;--fg:#e6edf3;--mut:#8b949e;--line:#30363d;--card:#161b22;--acc:#58a6ff;--ok:#3fb950;--okbg:#12261e;--ko:#d29922;--kobg:#272115}}
body{margin:0;background:var(--bg);color:var(--fg);font:16px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1040px;margin:0 auto;padding:32px 16px 64px}
h1{font-size:1.7rem;margin:0 0 4px} h2{margin-top:2.4rem;border-bottom:1px solid var(--line);padding-bottom:6px}
.sub{color:var(--mut);margin:0 0 20px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px 14px}
.kpi b{display:block;font-size:1.6rem} .kpi span{color:var(--mut);font-size:.85rem}
.scroll{overflow-x:auto} table{border-collapse:collapse;font-size:.88rem;margin:6px 0 14px}
th,td{padding:6px 10px;border:1px solid var(--line);text-align:left;vertical-align:top} thead th{background:var(--card)}
td.ok{background:var(--okbg);color:var(--ok);font-weight:600} td.ko{background:var(--kobg);color:var(--ko);font-weight:600}
td.na{color:var(--mut)}
.heat td{text-align:center;font-variant-numeric:tabular-nums}
.note{background:var(--card);border-left:3px solid var(--acc);padding:10px 14px;border-radius:0 6px 6px 0;margin:14px 0}
.legend{color:var(--mut);font-size:.85rem}
details summary{cursor:pointer;color:var(--acc)} img{max-width:100%;background:#fff;border-radius:6px}
"""

doc = f"""<!doctype html><html lang="fr"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Thymus signalisation stroma</title><style>{CSS}</style></head><body><main>
<h1>Thymus humain – signalisation stroma ↔ cellules T</h1>
<p class="sub">Reproduction du paragraphe « signalling factors known to be involved in thymic development » et de la fig. S11 de
Park et al., <i>Science</i> 367, eaay3224 (2020) · fichier HTA08_v01_A05_Science_human_fig1.h5ad (255 901 cellules) · CellPhoneDB 5.0.1</p>

<div class="kpis">
<div class="kpi"><b>{sum(v for *_, v, _ in CLAIMS)}/4</b><span>affirmations du papier retrouvées</span></div>
<div class="kpi"><b>99</b><span>groupes type × stade (≥50 cellules)</span></div>
<div class="kpi"><b>{len(ps_cols)}</b><span>couples testés (même stade)</span></div>
<div class="kpi"><b>{len(s5k & oursig)}/{len(s5k)}</b><span>paires de la table S5 significatives chez nous</span></div>
</div>

<h2>1. Reproduction de la fig. S11</h2>
<p>Comme dans le papier : types cellulaires (<code>Anno_level_fig1</code>) séparés par stade
(early 7–11 pcw, middle 12 pcw–10 mois, late ≥15 mois), groupes de moins de 50 cellules omis.
Toutes les cellules sont utilisées. Pourcentages = part des cellules du groupe qui expriment le gène.</p>
<div class="scroll"><table><thead><tr><th>Axe</th><th>Ce que dit le papier</th><th>Ce que l'on mesure</th><th>Verdict</th><th>Nuance</th></tr></thead>
<tbody>{claim_rows}</tbody></table></div>
<details><summary>Voir le dot plot (taille = % de cellules qui expriment, couleur = moyenne rapportée au maximum du gène)</summary>
<img alt="Dot plot des gènes de la fig. S11 par type cellulaire et stade" src="data:image/png;base64,{png}"></details>

<h2>2. CellPhoneDB par stade</h2>
<p>Même protocole que les auteurs (≤1000 cellules par groupe, 1000 permutations), avec un groupe par type × stade
et une analyse CellPhoneDB par stade (une analyse unique de 99 groupes dépassait la mémoire disponible).
Valeurs = p-value de CellPhoneDB, sender et receiver au même stade. Vert : p &lt; 0,05. « absent » : groupe de moins de 50 cellules.</p>
{step2}

<h2>3. Vue d'ensemble, stades regroupés</h2>
<p class="legend">Nombre de paires ligand–récepteur significatives (p &lt; 0,05) avec le ligand dans la ligne et le récepteur dans la colonne.
Seules les interactions de type ligand–récepteur sont comptées, orientées d'après la base (corrigé : la version précédente
mélangeait les paires d'adhésion et ne tenait pas compte des 77 paires où le récepteur est en premier).</p>
<div class="scroll">{heat}</div>

<h2>4. Comparaison avec la table S5</h2>
<ul>
<li>Comparaison par identifiants des partenaires, sans tenir compte de l'ordre (la table S5 écrit souvent
<code>IL7 receptor_IL7</code> là où notre base écrit <code>IL7_IL7_receptor</code>).</li>
<li>{len(s5k)} paires dans la table S5 ; {len(s5k & ours)} existent dans notre base v5 ; {len(s5k & oursig)} y sont significatives.</li>
<li>Les paires absentes viennent de l'ancienne version de CellPhoneDB (par exemple <code>CXCL12_CXCR3</code>, récepteurs ACKR).</li>
<li>Reproduction exacte impossible avec ce fichier : la table S5 couvre environ 53 types, dont certains absents du h5ad
(<code>NCC</code>, <code>Muscle_progenitors</code>, <code>Epi_PAX8</code>, <code>Fb_early</code>, <code>TEC_early</code>).</li>
</ul>

<h2>Limites et suite</h2>
<ul>
<li>Le stade est très lié au type cellulaire : cTEC 96 % early, mTEC(I) 79 % late, endothélium 92 % late. Les comparaisons entre stades ne valent que pour quelques types.</li>
<li>Petits groupes : ETP 316 cellules, endothélium lymphatique 115.</li>
<li>Moitié des cellules en chimie 3', moitié en 5' : à inclure comme variable de contexte dans MIIC.</li>
</ul>
<div class="note"><b>Étape suivante : CausalCCC, un axe à la fois.</b> 1) FGF7:FGFR2, fibroblastes → TEC ; 2) RANKL, ILC3 → mTEC(II) ;
3) Notch, TEC → ETP ; 4) LTB:LTBR en dernier (signal diffus). Pour MIIC, toutes les cellules des types choisis, avec le stade et la chimie comme variables de contexte.</div>
</main></body></html>"""

open(f'{RES}/first_results.html', 'w', encoding='utf-8').write(doc)
print('written', len(doc), 'bytes')
for a, _, c, v, _ in CLAIMS:
    print(a, 'OK' if v else 'CHECK', '|', c[:300])
