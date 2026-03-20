/**
 * docking.js - Épica 3: Molecular Docking & Drug Ranking
 * Lee scores reales de AutoDock Vina desde structural_metrics.json.
 * visualizePose carga complejos proteína+ligando cuando están disponibles.
 */

window.dockingViewer = null;
let isDockingInit = false;
let dockingUpdateLock = false;

window.initDockingView = async () => {
    if (!window.structureManifest) {
        try {
            const resp = await fetch(`data/manifest_structures.json?v=${Date.now()}`);
            window.structureManifest = await resp.json();
        } catch (e) { console.error("Error loading manifest:", e); }
    }

    const container = document.getElementById('docking-molstar-container');
    if (!isDockingInit && container) {
        isDockingInit = true;
        try {
            if (typeof PDBeMolstarPlugin !== 'undefined') {
                window.dockingViewer = new PDBeMolstarPlugin();
                await window.dockingViewer.render(container, {
                    customData: { url: 'data/structures/placeholder.pdb', format: 'pdb' },
                    bgColor: { r: 18, g: 18, b: 18 },
                    hideControls: true,
                    subscribeToEvents: false
                });
            }
        } catch (e) {
            console.error("Docking Viewer Init Error:", e);
            isDockingInit = false;
        }
    }

    loadDockingResults();
};

async function loadDockingResults() {
    const rankingContainer = document.getElementById('docking-ranking-list');
    const reproContainer   = document.getElementById('reposition-list');
    const scoreEl          = document.getElementById('top-vina-score');
    const reproCountEl     = document.getElementById('repro-count');
    const labelEl          = document.getElementById('docking-variant-label');

    if (!rankingContainer || !reproContainer) return;

    try {
        const resp = await fetch(`data/structural_metrics.json?v=${Date.now()}`);
        const metrics = await resp.json();

        // Seleccionar variante activa
        let variantId = window.selectedVariantId;
        if (!variantId || !metrics[variantId] || metrics[variantId].rmsd_active === '--') {
            const withData = Object.entries(metrics)
                .filter(([, m]) => m.rmsd_active !== '--' && !isNaN(parseFloat(m.rmsd_active)))
                .sort(([, a], [, b]) => Math.abs(parseFloat(b.rmsd_active)) - Math.abs(parseFloat(a.rmsd_active)));
            variantId = withData.length > 0 ? withData[0][0] : Object.keys(metrics)[0];
        }

        const parts  = variantId.split('_');
        const hugo   = parts[0];
        const change = parts.slice(1).join('_');
        if (labelEl) labelEl.textContent = `${hugo} ${change || 'WT'}`;

        const m = metrics[variantId] || {};

        // Determinar fuente de scores
        const hasVinaScores = m.vina_scores && Object.keys(m.vina_scores).length > 0;
        const rmsdActive    = parseFloat(m.rmsd_active);
        const hasRmsdData   = !isNaN(rmsdActive);

        const drugCatalog = [
            { name: "Alpelisib",    target: "PI3Kα",     status: "FDA",   source: "ChEMBL",   baseWT: -8.2, selectivity: "high"   },
            { name: "Abemaciclib",  target: "CDK4/6",    status: "FDA",   source: "DrugBank",  baseWT: -8.3, selectivity: "high"   },
            { name: "Palbociclib",  target: "CDK4/6",    status: "FDA",   source: "DrugBank",  baseWT: -8.1, selectivity: "high"   },
            { name: "Ribociclib",   target: "CDK4/6",    status: "FDA",   source: "DrugBank",  baseWT: -8.0, selectivity: "high"   },
            { name: "Copanlisib",   target: "Pan-PI3K",  status: "FDA",   source: "ChEMBL",   baseWT: -7.9, selectivity: "medium" },
            { name: "Capivasertib", target: "AKT1/2/3", status: "FDA",   source: "DrugBank",  baseWT: -7.4, selectivity: "high"   },
            { name: "Everolimus",   target: "mTOR",      status: "FDA",   source: "ChEMBL",   baseWT: -7.6, selectivity: "medium" },
            { name: "Gedatolisib",  target: "PI3K/mTOR", status: "Trial", source: "ChEMBL",   baseWT: -7.8, selectivity: "medium" },
            { name: "Taselisib",    target: "PI3Kα",     status: "Trial", source: "ChEMBL",   baseWT: -7.5, selectivity: "high"   },
            { name: "Pictilisib",   target: "PI3K",      status: "Trial", source: "ChEMBL",   baseWT: -7.2, selectivity: "medium" },
            { name: "Buparlisib",   target: "Pan-PI3K",  status: "Trial", source: "DrugBank",  baseWT: -6.9, selectivity: "low"    },
        ];

        let results, dataSource;

        if (hasVinaScores) {
            dataSource = "vina_real";
            results = drugCatalog.map(d => {
                const score   = m.vina_scores[d.name];
                const affinity = score !== undefined ? score : d.baseWT;
                const delta   = parseFloat((affinity - d.baseWT).toFixed(1));
                return { ...d, affinity, delta, isReal: score !== undefined };
            }).sort((a, b) => a.affinity - b.affinity);
        } else {
            dataSource = hasRmsdData ? "rmsd_calc" : "baseline";
            const impactFactor = hasRmsdData ? (Math.abs(rmsdActive) - 2.0) * 0.3 : 0;
            const mult = { high: 1.2, medium: 0.8, low: 0.4 };
            results = drugCatalog.map(d => {
                const adj   = parseFloat((d.baseWT - impactFactor * (mult[d.selectivity] || 1)).toFixed(1));
                const delta = parseFloat((adj - d.baseWT).toFixed(1));
                return { ...d, affinity: adj, delta, isReal: false };
            }).sort((a, b) => a.affinity - b.affinity);
        }

        // Verificar qué fármacos tienen complejo disponible
        const complexAvailable = {};
        await Promise.all(results.map(async d => {
            const drugSafe = d.name.replace(/ /g, '_');
            const url = `data/docking/complex/${variantId}_${drugSafe}_complex.pdb`;
            try {
                const r = await fetch(url, { method: 'HEAD' });
                complexAvailable[d.name] = r.ok;
            } catch { complexAvailable[d.name] = false; }
        }));

        // Render ranking
        const statusClass = { FDA: 'ds-fda', Trial: 'ds-trial', Exp: 'ds-exp' };
        let rankingHtml = '';

        results.slice(0, 10).forEach((drug, idx) => {
            const deltaStr   = drug.delta < 0 ? `${drug.delta}` : drug.delta > 0 ? `+${drug.delta}` : '—';
            const deltaColor = drug.delta < -0.3 ? 'var(--success)' : drug.delta > 0.3 ? 'var(--danger)' : 'var(--text-tertiary)';
            const realBadge  = drug.isReal
                ? `<span style="font-size:0.6rem;background:rgba(16,185,129,0.15);color:var(--success);padding:1px 5px;border-radius:3px;margin-left:4px;">Vina</span>`
                : '';
            const hasComplex = complexAvailable[drug.name];
            const poseTitle  = hasComplex ? '🔬 Pose (proteína+ligando)' : '🔬 Pose (proteína)';
            const poseBg     = hasComplex
                ? 'background:rgba(13,115,119,0.25);border-color:var(--accent);color:var(--accent);'
                : 'background:rgba(13,115,119,0.1);border-color:var(--accent);color:var(--accent);';

            rankingHtml += `
            <div style="display:flex;align-items:center;gap:10px;padding:10px 0;border-bottom:1px solid var(--border);flex-wrap:nowrap;">
                <div class="drug-rank ${idx < 3 ? 'top' : ''}" style="flex-shrink:0;">${idx + 1}</div>
                <div style="flex:1;min-width:0;">
                    <div style="font-weight:600;font-size:0.85rem;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;">
                        ${drug.name}${realBadge}
                    </div>
                    <div style="font-size:0.7rem;color:var(--text-tertiary);margin-top:2px;">
                        ${drug.target} · <span style="opacity:0.6;">${drug.source}</span>
                    </div>
                </div>
                <div style="font-family:monospace;font-size:0.82rem;font-weight:600;flex-shrink:0;white-space:nowrap;">
                    ${drug.affinity} kcal/mol
                </div>
                <div class="${statusClass[drug.status] || 'ds-exp'}" style="flex-shrink:0;">${drug.status}</div>
                <div style="font-size:0.8rem;font-weight:700;color:${deltaColor};flex-shrink:0;min-width:36px;text-align:right;">
                    ${deltaStr}
                </div>
                <button
                    onclick="visualizePose('${variantId}', '${drug.name}')"
                    title="${poseTitle}"
                    style="flex-shrink:0;padding:5px 10px;font-size:0.7rem;${poseBg}border-radius:6px;border:1px solid;cursor:pointer;white-space:nowrap;transition:all 0.2s;"
                    onmouseover="this.style.background='var(--accent)';this.style.color='#fff';"
                    onmouseout="this.style.background='${hasComplex ? 'rgba(13,115,119,0.25)' : 'rgba(13,115,119,0.1)'}';this.style.color='var(--accent)';"
                >🔬 Pose${hasComplex ? ' ⚛' : ''}</button>
            </div>`;
        });

        rankingContainer.innerHTML = rankingHtml || '<div style="padding:20px;color:var(--text-tertiary);">Sin datos.</div>';

        // Render repositionables
        const repositionables = results.filter(d => d.delta <= -0.5);
        let reproHtml = '';
        repositionables.forEach(drug => {
            const hasC = complexAvailable[drug.name];
            reproHtml += `
            <div class="reposition-card">
                <div style="font-size:0.85rem;font-weight:600;margin-bottom:8px;">
                    ${drug.name} <span style="color:var(--text-tertiary);font-weight:400;">— ${drug.target}</span>
                    ${drug.isReal ? '<span style="font-size:0.65rem;background:rgba(16,185,129,0.15);color:var(--success);padding:1px 5px;border-radius:3px;margin-left:4px;">Vina</span>' : ''}
                </div>
                <div style="display:flex;justify-content:space-between;font-size:0.75rem;margin-bottom:4px;">
                    <span style="color:var(--text-tertiary);">Baseline WT</span>
                    <span style="font-weight:600;">${drug.baseWT} kcal/mol</span>
                </div>
                <div style="display:flex;justify-content:space-between;font-size:0.75rem;margin-bottom:4px;">
                    <span style="color:var(--text-tertiary);">Score ${change}</span>
                    <span style="font-weight:600;color:var(--success);">${drug.affinity} kcal/mol</span>
                </div>
                <div style="display:flex;justify-content:space-between;font-size:0.75rem;margin-bottom:10px;">
                    <span style="color:var(--text-tertiary);">Δ Impacto</span>
                    <span class="delta-neg">${drug.delta} kcal/mol</span>
                </div>
                <button onclick="visualizePose('${variantId}', '${drug.name}')"
                    class="fasta-btn active"
                    style="width:100%;font-size:0.75rem;padding:6px;">
                    🔬 Ver Pose${hasC ? ' (proteína+ligando) ⚛' : ''}
                </button>
                <div class="reposition-badge">✓ Candidato repositionable</div>
            </div>`;
        });

        reproContainer.innerHTML = reproHtml ||
            `<div style="padding:20px;font-size:0.8rem;color:var(--text-tertiary);text-align:center;">
                Sin candidatos repositionables para <strong>${hugo} ${change}</strong>.
            </div>`;

        if (scoreEl)      scoreEl.textContent     = results[0]?.affinity ?? '--';
        if (reproCountEl) reproCountEl.textContent = repositionables.length;

        // Nota de fuente
        const noteId = 'docking-data-note';
        let note = document.getElementById(noteId);
        if (!note && rankingContainer.parentElement) {
            note = document.createElement('div');
            note.id = noteId;
            note.style.cssText = 'font-size:0.7rem;color:var(--text-tertiary);margin-top:10px;padding-top:8px;border-top:1px solid var(--border);';
            rankingContainer.parentElement.appendChild(note);
        }
        if (note) {
            const complexCount = Object.values(complexAvailable).filter(Boolean).length;
            const srcMsg = {
                vina_real: `✓ Scores reales AutoDock Vina 1.2.3 · Badge "Vina" = score medido`,
                rmsd_calc: `⚙ Scores estimados desde RMSD (${rmsdActive?.toFixed(2)} Å)`,
                baseline:  `ℹ Scores de referencia — sin datos estructurales`,
            };
            note.textContent = `${srcMsg[dataSource] || ''} · ⚛ = complejo proteína+ligando disponible (${complexCount}/${results.length})`;
        }

    } catch (err) {
        console.error("Error en loadDockingResults:", err);
        if (rankingContainer) rankingContainer.innerHTML =
            `<div style="color:var(--danger);padding:20px;">Error: ${err.message}</div>`;
    }
}

function safeSetText(id, text) {
    const el = document.getElementById(id);
    if (el) el.textContent = text;
}

window.visualizePose = async (variantId, drugName) => {
    if (!window.dockingViewer?.visual || !window.structureManifest || dockingUpdateLock) return;

    const drugSafe = drugName.replace(/ /g, '_');

    // Prioridad de carga:
    // 1. Complejo proteína+ligando (mejor vista)
    // 2. Solo ligando _out.pdbqt (pose Vina)
    // 3. Solo receptor (estructura proteína)
    const complexUrl  = `data/docking/complex/${variantId}_${drugSafe}_complex.pdb`;
    const ligandUrl   = `data/docking/results/${variantId}_${drugSafe}_out.pdbqt`;
    const receptorUrl = window.structureManifest[variantId]
                     || window.structureManifest[variantId.split('_')[0] + '_WT'];

    let urlToLoad  = receptorUrl || null;
    let viewMode   = 'receptor';

    try {
        const r1 = await fetch(complexUrl, { method: 'HEAD' });
        if (r1.ok) { urlToLoad = complexUrl; viewMode = 'complex'; }
        else {
            const r2 = await fetch(ligandUrl, { method: 'HEAD' });
            if (r2.ok) { urlToLoad = ligandUrl; viewMode = 'ligand'; }
        }
    } catch (e) { /* usar fallback */ }

    if (!urlToLoad) {
        alert(`Sin estructura disponible para ${variantId}.`);
        return;
    }

    // Scroll al viewer
    document.getElementById('docking-viewer-container')
        ?.scrollIntoView({ behavior: 'smooth', block: 'center' });

    // Indicador de carga
    const sub = document.querySelector('#view-docking .view-subtitle');
    if (sub) {
        const modeLabel = {
            complex:  `⚛ Complejo proteína+ligando`,
            ligand:   `💊 Pose del ligando (Vina pose 1)`,
            receptor: `🧬 Estructura proteica (sin ligando)`
        };
        sub.textContent = `${modeLabel[viewMode]}: ${drugName} → ${variantId.replace('_', ' ')} | AutoDock Vina 1.2.3`;
    }

    dockingUpdateLock = true;
    try {
        const format = urlToLoad.endsWith('.cif') ? 'mmcif' : 'pdb';
        await window.dockingViewer.visual.update({ customData: { url: urlToLoad, format } });

        // En modo complejo: resaltar ligando (cadena L) y sitio de mutación
        if (viewMode === 'complex') {
            const mutMatch = variantId.match(/[A-Za-z](\d+)[A-Za-z]/);
            setTimeout(async () => {
                try {
                    const selectData = [
                        // Ligando en naranja brillante
                        { chain_id: 'L', color: { r: 255, g: 140, b: 0 }, focus: true }
                    ];
                    // Residuo mutado en rojo
                    if (mutMatch) {
                        selectData.push({
                            chain_id: 'A',
                            residue_number: parseInt(mutMatch[1]),
                            color: { r: 239, g: 68, b: 68 }
                        });
                    }
                    await window.dockingViewer.visual.select({ data: selectData });
                } catch (e) { /* silencioso */ }
            }, 700);
        } else if (viewMode === 'receptor') {
            // Solo resaltar residuo mutado
            const mutMatch = variantId.match(/[A-Za-z](\d+)[A-Za-z]/);
            if (mutMatch) {
                setTimeout(async () => {
                    try {
                        await window.dockingViewer.visual.select({
                            data: [{ residue_number: parseInt(mutMatch[1]), color: { r: 255, g: 165, b: 0 }, focus: true }]
                        });
                    } catch (e) {}
                }, 600);
            }
        }

    } catch (e) {
        console.warn("Error loading pose:", e);
    } finally {
        setTimeout(() => { dockingUpdateLock = false; }, 100);
    }
};
