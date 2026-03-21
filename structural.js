/**
 * structural.js - Análisis Estructural 3D Comparativo
 * PDBe Molstar v3.2.0
 *
 * ESTRATEGIA FINAL:
 * render() con customData usando un CIF local real (CBFB_WT, 182 aa).
 * Este archivo existe siempre en el proyecto y es mmCIF válido.
 * visual.update() lo reemplaza con la estructura seleccionada.
 */

window.wtViewer  = null;
window.mutViewer = null;
let structuralMetrics   = {};
let isGlobalViewerInit  = false;
let isScientificMode    = false;
let linkedViews         = true;
let isSuperpositionMode = false;
let viewerUpdateLock    = { wt: false, mut: false };
let _activeVariant      = null;

// ── Silenciar warnings internos ───────────────────────────────────────────
const _origWarn  = console.warn;
const _origError = console.error;
const _skipMsgs  = ['already added. Call removeSymbol', 'Invalid PDBeMolstarPlugin'];
console.warn  = (...a) => { if (typeof a[0]==='string' && _skipMsgs.some(s=>a[0].includes(s))) return; _origWarn.apply(console,a); };
console.error = (...a) => { if (typeof a[0]==='string' && _skipMsgs.some(s=>a[0].includes(s))) return; _origError.apply(console,a); };

// ── Inicialización ────────────────────────────────────────────────────────
window.initStructuralView = async () => {

    if (Object.keys(structuralMetrics).length === 0) {
        try {
            const [mr, mfr] = await Promise.all([
                fetch(`data/structural_metrics.json?v=${Date.now()}`),
                fetch(`data/manifest_structures.json?v=${Date.now()}`)
            ]);
            structuralMetrics        = await mr.json();
            window.structureManifest = await mfr.json();
        } catch(e) { console.error('Error loading data:', e); }
    }

    const wtEl  = document.getElementById('wt-molstar-container');
    const mutEl = document.getElementById('mut-molstar-container');

    if (!isGlobalViewerInit && wtEl && mutEl) {
        isGlobalViewerInit = true;
        try {
            // Usar un CIF local real y válido para el render() inicial.
            // CBFB_WT.cif es un Boltz mmCIF de 182 aa — siempre disponible.
            const initUrl = 'data/structures/boltz_results_CBFB_WT/predictions/CBFB_WT/CBFB_WT_model_0.cif';
            const cfg = {
                customData:   { url: initUrl, format: 'mmcif' },
                bgColor:      { r: 18, g: 18, b: 18 },
                hideControls: true,
                landscape:    false,
            };

            window.wtViewer  = new PDBeMolstarPlugin();
            window.mutViewer = new PDBeMolstarPlugin();

            // Renderizar secuencialmente para evitar conflictos WebGL
            await window.wtViewer.render(wtEl,  cfg);
            await window.mutViewer.render(mutEl, cfg);

            const ro = new ResizeObserver(() => window.dispatchEvent(new Event('resize')));
            ro.observe(wtEl);
            ro.observe(mutEl);

            setTimeout(setupLinkedViews, 1500);

        } catch(e) {
            console.error('Molstar init error:', e);
            isGlobalViewerInit = false;
        }
    }

    // Botón Modo Científico
    const sciBtn = document.getElementById('btn-scientific-mode');
    if (sciBtn && !sciBtn._bound) {
        sciBtn._bound = true;
        sciBtn.onclick = () => {
            isScientificMode = !isScientificMode;
            sciBtn.textContent = `🔬 Modo Científico: ${isScientificMode ? 'ON' : 'OFF'}`;
            sciBtn.className   = isScientificMode ? 'fasta-btn active' : 'fasta-btn';
            const sel = document.getElementById('struct-variant-select');
            if (sel?.value) updateStructuralDisplay(sel.value);
        };
    }

    // Selector de variantes
    const select = document.getElementById('struct-variant-select');
    if (select && select.options.length <= 1) {
        try {
            const r = await fetch(`data/top_variants.json?v=${Date.now()}`);
            const variants = await r.json();
            select.innerHTML = '<option value="">-- Cambiar Variante --</option>';
            variants.forEach(v => {
                const o = document.createElement('option');
                o.value       = `${v.hugo_symbol}_${v.protein_change}`;
                o.textContent = `${v.hugo_symbol} ${v.protein_change}`;
                select.appendChild(o);
            });
        } catch(e) {}
    }

    if (select) {
        if (window.selectedVariantId) select.value = window.selectedVariantId;
        else if (select.options.length > 1) select.selectedIndex = 1;

        // Esperar a que los viewers estén listos antes de cargar la estructura real
        waitForBothViewers(8000).then(() => {
            if (select.value) updateStructuralDisplay(select.value);
        });

        select.onchange = e => {
            _activeVariant = e.target.value;
            updateStructuralDisplay(e.target.value);
        };
    }
};

// ── Esperar a que ambos viewers tengan canvas3d ───────────────────────────
function waitForBothViewers(timeout = 8000) {
    return new Promise(resolve => {
        const t0 = Date.now();
        const id = setInterval(() => {
            if (window.wtViewer?.plugin?.canvas3d && window.mutViewer?.plugin?.canvas3d) {
                clearInterval(id); resolve();
            } else if (Date.now() - t0 > timeout) {
                clearInterval(id); resolve();
            }
        }, 200);
    });
}

const waitForPlugin = (plugin, timeout = 8000) => new Promise((resolve, reject) => {
    if (plugin?.plugin?.canvas3d) { resolve(); return; }
    const t0 = Date.now();
    const id = setInterval(() => {
        if (plugin?.plugin?.canvas3d)       { clearInterval(id); resolve(); }
        else if (Date.now() - t0 > timeout) { clearInterval(id); reject(new Error('timeout')); }
    }, 150);
});

// ── Actualizar display ────────────────────────────────────────────────────
async function updateStructuralDisplay(variantId) {
    if (!variantId) return;
    _activeVariant = variantId;

    const parts  = variantId.split('_');
    const hugo   = parts[0];
    const change = parts.slice(1).join('_');

    safeSetText('struct-view-title', `Análisis estructural — ${hugo} ${change}`);
    safeSetText('mut-panel-title',   `Proteína mutada (${change})`);

    const m = structuralMetrics[variantId] || {};
    safeSetText('wt-avg-plddt',  m.wt_plddt  || '--');
    safeSetText('rmsd-global',   m.rmsd_global && m.rmsd_global !== '--' ? m.rmsd_global + ' Å' : '--');
    safeSetText('rmsd-active',   m.rmsd_active && m.rmsd_active !== '--' ? m.rmsd_active + ' Å' : '--');
    safeSetText('wt-confidence', parseFloat(m.wt_plddt) > 70 ? 'Alta' : parseFloat(m.wt_plddt) > 50 ? 'Media' : '--');

    const wtBar  = document.getElementById('wt-plddt-bar');
    const mutBar = document.getElementById('mut-plddt-bar');
    if (wtBar)  wtBar.style.width = `${Math.min(parseFloat(m.wt_plddt)  || 0, 100)}%`;
    if (mutBar) mutBar.style.width = `${Math.min(parseFloat(m.mut_plddt) || 0, 100)}%`;

    if (m.rmsd_regions?.length) renderRMSDChart(m.rmsd_regions);

    try {
        await Promise.all([
            waitForPlugin(window.wtViewer),
            waitForPlugin(window.mutViewer)
        ]);
        if (_activeVariant !== variantId) return;
        await loadBestStructure(window.wtViewer,  `${hugo}_WT`, hugo, 'wt',  variantId);
        await loadBestStructure(window.mutViewer, variantId,    hugo, 'mut', variantId);
    } catch(e) {
        console.warn('[Molstar] Reintentando en 2s');
        setTimeout(() => { if (_activeVariant === variantId) updateStructuralDisplay(variantId); }, 2000);
    }
}

// ── Cargar estructura ─────────────────────────────────────────────────────
async function loadBestStructure(plugin, id, hugo, lockKey, requestedVariant) {
    if (!plugin?.visual || !window.structureManifest || viewerUpdateLock[lockKey]) return;
    if (_activeVariant !== requestedVariant) return;

    document.querySelector(`#${lockKey}-viewer-panel .no-structure-overlay`)?.remove();

    const url = window.structureManifest[id];
    if (!url) { setViewerUnavailable(lockKey, id); return; }

    viewerUpdateLock[lockKey] = true;
    setViewerLoading(lockKey, true);
    try {
        const format = url.endsWith('.cif') ? 'mmcif' : 'pdb';
        await plugin.visual.update({ customData: { url: `${url}?v=${Date.now()}`, format } });

        if (_activeVariant !== requestedVariant) return;
        if (isScientificMode) applyScientificFidelity(plugin);

        if (lockKey === 'mut') {
            const match = id.match(/[A-Z](\d+)[A-Z]/);
            if (match) {
                const resNum = parseInt(match[1]);
                setTimeout(() => {
                    if (_activeVariant === requestedVariant)
                        highlightMutationSite(plugin, resNum);
                }, 600);
            }
        }
    } catch(e) {
        console.warn(`[Molstar] Load issue ${id}:`, e);
        setViewerUnavailable(lockKey, id);
    } finally {
        setViewerLoading(lockKey, false);
        setTimeout(() => { viewerUpdateLock[lockKey] = false; }, 150);
    }
}

// ── UI helpers ────────────────────────────────────────────────────────────
function setViewerLoading(key, loading) {
    const el = document.getElementById(`${key}-plddt-label`);
    if (!el) return;
    if (loading) { el.textContent = '⏳ Cargando...'; el.style.color = 'var(--accent)'; }
    else {
        const p = document.getElementById(`${key}-avg-plddt`)?.textContent;
        el.textContent = `pLDDT: ${p || '--'}`;
        el.style.color = 'rgba(255,255,255,0.7)';
    }
}

function setViewerUnavailable(key, id) {
    const el = document.getElementById(`${key}-plddt-label`);
    if (el) { el.textContent = '⚠ No disponible'; el.style.color = '#f87171'; }
    const panel = document.getElementById(`${key}-viewer-panel`);
    if (panel && !panel.querySelector('.no-structure-overlay')) {
        const ov = document.createElement('div');
        ov.className = 'no-structure-overlay';
        ov.innerHTML = `<div style="text-align:center;color:rgba(255,255,255,0.5);">
            <div style="font-size:2.5rem;margin-bottom:12px;">🧬</div>
            <div style="font-size:0.9rem;font-weight:600;color:#fff;">Estructura no disponible</div>
            <div style="font-size:0.75rem;margin-top:6px;font-family:monospace;opacity:0.7;">${id}</div>
        </div>`;
        panel.appendChild(ov);
    }
}

// ── Representaciones ──────────────────────────────────────────────────────
const REPR_MAP = {
    cartoon:   { type: 'cartoon',           color: 'sequence-id' },
    surface:   { type: 'molecular-surface', color: 'hydrophobicity' },
    ballstick: { type: 'ball-and-stick',    color: 'element-symbol' },
};
window.setRepresentation = (key, reprType) => {
    const p = { wt: window.wtViewer, mut: window.mutViewer }[key];
    if (!p?.visual) return;
    p.visual.update({ visual: REPR_MAP[reprType] });
    document.querySelector(`#${key}-viewer-panel .viewer-toolbar`)
        ?.querySelectorAll('.t-btn')
        .forEach(b => b.classList.toggle('active', b.textContent.toLowerCase().includes(reprType.slice(0,4))));
};

function applyScientificFidelity(plugin) {
    if (!plugin?.visual) return;
    plugin.visual.update({ visual: { color: 'plddt-confidence' } })
        .catch(() => plugin.visual.update({ visual: { color: 'b-factor' } }));
}

async function highlightMutationSite(plugin, residueNumber) {
    if (!plugin?.visual) return;
    try { await plugin.visual.select({ data: [{ residue_number: residueNumber, color: { r:239, g:68, b:68 }, focus: true }] }); }
    catch(e) {}
}

// ── Superposición ─────────────────────────────────────────────────────────
window.toggleSuperposition = async () => {
    const sel = document.getElementById('struct-variant-select');
    const vid = sel?.value;
    const btn = document.getElementById('btn-superposition');
    if (!vid) {
        if (btn) {
            btn.style.borderColor = 'var(--danger)';
            btn.style.color = 'var(--danger)';
            setTimeout(() => {
                btn.style.borderColor = 'var(--info)';
                btn.style.color = 'var(--info)';
            }, 1500);
        }
        return;
    }
    const hugo = vid.split('_')[0];
    if (isSuperpositionMode) {
        isSuperpositionMode = false;
        if (btn) {
            btn.style.background = '';
            btn.style.color = 'var(--info)';
            btn.style.borderColor = 'var(--info)';
            btn.textContent = '🔄 Superponer WT/Mut';
        }
        safeSetText('mut-panel-title', 'Proteína mutada');
        await loadBestStructure(window.mutViewer, vid, hugo, 'mut', vid);
        return;
    }
    const wtUrl = window.structureManifest?.[`${hugo}_WT`];
    if (!wtUrl) return;
    isSuperpositionMode = true;
    if (btn) {
        btn.style.background = 'var(--info)';
        btn.style.color = '#fff';
        btn.textContent = '🔄 Superposición: ON';
    }
    safeSetText('mut-panel-title', 'Superposición: WT');
    setViewerLoading('mut', true);
    try {
        await window.mutViewer.visual.update({
            customData: { url: wtUrl, format: wtUrl.endsWith('.cif') ? 'mmcif' : 'pdb' }
        });
    } catch(e) {
        isSuperpositionMode = false;
        if (btn) {
            btn.style.background = '';
            btn.style.color = 'var(--info)';
            btn.textContent = '🔄 Superponer WT/Mut';
        }
    } finally {
        setViewerLoading('mut', false);
    }
};

// ── Linked views ──────────────────────────────────────────────────────────
function setupLinkedViews(retries = 15) {
    const wm = window.wtViewer?.plugin;
    const mm = window.mutViewer?.plugin;
    if (!wm?.canvas3d || !mm?.canvas3d) {
        if (retries > 0) setTimeout(() => setupLinkedViews(retries - 1), 500);
        return;
    }
    wm.canvas3d.camera.stateChanged.subscribe(s => { if (linkedViews) mm.canvas3d.camera.setState(s, 0); });
    ['wt-sync-status','mut-sync-status'].forEach(sid => {
        const el = document.getElementById(sid);
        if (!el) return;
        el.textContent = '🔗 SinC: ON'; el.style.cursor = 'pointer';
        el.onclick = () => {
            linkedViews = !linkedViews;
            ['wt-sync-status','mut-sync-status'].forEach(s => {
                const e = document.getElementById(s);
                if (e) { e.textContent = linkedViews ? '🔗 SinC: ON' : '🔗 SinC: OFF'; e.classList.toggle('off', !linkedViews); }
            });
        };
    });
}

// ── Export PNG ────────────────────────────────────────────────────────────
window.exportViewerImage = key => {
    const canvas = document.getElementById(`${key}-molstar-container`)?.querySelector('canvas');
    if (!canvas) return;
    const a = document.createElement('a');
    a.download = `${document.getElementById('struct-variant-select')?.value || 'structure'}_${key}.png`;
    a.href = canvas.toDataURL('image/png');
    a.click();
};

// ── RMSD Chart ────────────────────────────────────────────────────────────
function renderRMSDChart(regions) {
    const c = document.getElementById('rmsd-chart-container');
    if (!c) return;
    c.innerHTML = '';
    regions.forEach(reg => {
        const row = document.createElement('div');
        row.style.cssText = 'display:flex;align-items:center;gap:15px;margin-bottom:8px;';
        row.innerHTML = `
            <div style="width:150px;font-size:0.85rem;font-weight:500;color:var(--text-secondary);">${reg.region}</div>
            <div style="flex:1;height:10px;background:rgba(255,255,255,0.05);border-radius:5px;overflow:hidden;">
                <div style="width:${Math.min((Math.abs(reg.value)/10)*100,100)}%;height:100%;background:${reg.color||'var(--accent)'};border-radius:5px;"></div>
            </div>
            <div style="width:50px;font-size:0.8rem;text-align:right;color:var(--text-primary);">${reg.value} Å</div>`;
        c.appendChild(row);
    });
}

function safeSetText(id, text) {
    const el = document.getElementById(id);
    if (el) el.textContent = text;
}
