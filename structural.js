/**
 * structural.js - Stable Comparative Protein Structure Analysis
 */

let wtViewer = null;
let mutViewer = null;
let structuralMetrics = {};
let isGlobalViewerInit = false;

window.initStructuralView = async () => {
    // 1. Data Initialization
    if (Object.keys(structuralMetrics).length === 0) {
        try {
            const resp = await fetch('data/structural_metrics.json');
            structuralMetrics = await resp.json();
        } catch (e) { console.error("Error loading metrics:", e); }
    }

    // 2. Strict Singleton Initialization
    const wtContainer = document.getElementById('wt-molstar-container');
    const mutContainer = document.getElementById('mut-molstar-container');

    if (!isGlobalViewerInit) {
        // Set lock immediately to prevent race conditions during async render
        isGlobalViewerInit = true; 

        if (wtContainer && !wtViewer) {
            try {
                wtViewer = new PDBeMolstarPlugin();
                wtViewer.render(wtContainer, {
                    // CRITICAL: Initialize with a valid LOCAL placeholder to avoid EBI 404s
                    customData: { url: 'data/structures/placeholder.pdb', format: 'pdb' },
                    bgColor: { r: 18, g: 18, b: 18 },
                    hideControls: true,
                    subscribeToEvents: false
                });
            } catch (e) { console.error("WT Init Error:", e); }
        }
        
        // Slight delay between instances
        await new Promise(r => setTimeout(r, 600));

        if (mutContainer && !mutViewer) {
            try {
                mutViewer = new PDBeMolstarPlugin();
                mutViewer.render(mutContainer, {
                    customData: { url: 'data/structures/placeholder.pdb', format: 'pdb' },
                    bgColor: { r: 18, g: 18, b: 18 },
                    hideControls: true,
                    subscribeToEvents: false
                });
            } catch (e) { console.error("MUT Init Error:", e); }
        }
    }

    // 3. Populate Selector
    const select = document.getElementById('struct-variant-select');
    if (select && select.options.length <= 1) {
        try {
            const response = await fetch('data/top_variants.json');
            const variants = await response.json();
            
            select.innerHTML = '<option value="">-- Cambiar Variante --</option>';
            variants.forEach(v => {
                const opt = document.createElement('option');
                opt.value = `${v.hugo_symbol}_${v.protein_change}`;
                opt.textContent = `${v.hugo_symbol} ${v.protein_change}`;
                select.appendChild(opt);
            });
        } catch (err) { console.error("Error populating structural selector:", err); }
    }

    // Sync state
    if (select) {
        if (window.selectedVariantId) select.value = window.selectedVariantId;
        else if (select.options.length > 1) select.selectedIndex = 1;
        
        updateStructuralDisplay(select.value);
        select.onchange = (e) => updateStructuralDisplay(e.target.value);
    }
};

async function updateStructuralDisplay(variantId) {
    if (!variantId) return;
    const [hugo, change] = variantId.split('_');
    
    safeSetText('struct-view-title', `Análisis estructural — ${hugo} ${change}`);
    safeSetText('mut-panel-title', `Proteína mutada (${change})`);

    const metrics = structuralMetrics[variantId] || (structuralMetrics["PIK3CA_H1047R"] || {});

    safeSetText('wt-avg-plddt', metrics.wt_plddt || '--');
    safeSetText('wt-active-plddt', metrics.wt_plddt_active || '--');
    safeSetText('wt-plddt-label', `pLDDT: ${metrics.wt_plddt || '--'}`);
    
    if (document.getElementById('wt-plddt-bar')) {
        document.getElementById('wt-plddt-bar').style.width = `${metrics.wt_plddt || 0}%`;
    }
    
    safeSetText('wt-confidence', (metrics.wt_plddt > 80 ? 'Alta' : 'Media'));

    safeSetText('mut-avg-plddt', metrics.mut_plddt || '--');
    safeSetText('rmsd-global', (metrics.rmsd_global || '--') + ' Å');
    safeSetText('rmsd-active', (metrics.rmsd_active || '--') + ' Å');
    safeSetText('mut-plddt-label', `pLDDT: ${metrics.mut_plddt || '--'}`);
    
    if (document.getElementById('mut-plddt-bar')) {
        document.getElementById('mut-plddt-bar').style.width = `${metrics.mut_plddt || 0}%`;
    }

    if (metrics.rmsd_regions) renderRMSDChart(metrics.rmsd_regions);

    // 4. Load Models using DIRECT DATA to avoid fetch issues in plugin
    const wtUrl = `data/structures/${hugo}_WT.pdb`;
    const mutUrl = `data/structures/${hugo}_${change}_MUT.pdb`;
    
    loadModelRobustly(wtViewer, wtUrl);
    loadModelRobustly(mutViewer, mutUrl);
}

function safeSetText(id, text) {
    const el = document.getElementById(id);
    if (el) el.textContent = text;
}

function renderRMSDChart(regions) {
    const container = document.getElementById('rmsd-chart-container');
    if (!container) return;
    container.innerHTML = '';
    
    regions.forEach(reg => {
        const row = document.createElement('div');
        row.style.display = 'flex'; row.style.alignItems = 'center'; row.style.gap = '15px';
        row.innerHTML = `
            <div style="width: 150px; font-size: 0.85rem; font-weight: 500; color: var(--text-secondary); text-align: left;">${reg.region}</div>
            <div style="flex: 1; height: 10px; background: rgba(255,255,255,0.05); border-radius: 5px; overflow: hidden;">
                <div style="width: ${(reg.value / 5) * 100}%; height: 100%; background: ${reg.color};"></div>
            </div>
            <div style="width: 50px; font-size: 0.8rem; text-align: right; color: var(--text-primary);">${reg.value} Å</div>
        `;
        container.appendChild(row);
    });
}

async function loadModelRobustly(plugin, url) {
    if (!plugin || !plugin.visual) return;
    try {
        const resp = await fetch(url);
        if (!resp.ok) { console.warn("PDB mirror missing:", url); return; }
        
        const pdbData = await resp.text();
        
        // Exact pattern from working version: customData with literal 'data'
        await plugin.visual.update({
            customData: {
                data: pdbData,
                format: 'pdb'
            }
        }, true);
    } catch (e) {
        console.warn("Molstar silent error:", e.message);
    }
}

// Global upload handler
const uploadEl = document.getElementById('pdb-upload');
if (uploadEl) {
    uploadEl.onchange = (e) => {
        const file = e.target.files[0];
        if (!file || !mutViewer) return;
        const reader = new FileReader();
        reader.onload = async (evt) => {
            try {
                await mutViewer.visual.update({
                    customData: { data: evt.target.result, format: 'pdb' }
                }, true);
            } catch (err) { alert("Error cargando PDB."); }
        };
        reader.readAsText(file);
    };
}
