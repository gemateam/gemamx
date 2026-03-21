/**
 * dashboard.js — Dashboard de resumen ejecutivo
 * Lee top_variants.json, structural_metrics.json y manifest_structures.json
 * para mostrar KPIs reales del estado del pipeline.
 */

window.initDashboard = async () => {
    console.log("Iniciando Dashboard...");
    try {
        const [varResp, metResp, manResp] = await Promise.all([
            fetch(`data/top_variants.json?v=${Date.now()}`),
            fetch(`data/structural_metrics.json?v=${Date.now()}`),
            fetch(`data/manifest_structures.json?v=${Date.now()}`)
        ]);

        if (!varResp.ok || !metResp.ok || !manResp.ok) {
            throw new Error(`Error en fetch: ${varResp.status} ${metResp.status} ${manResp.status}`);
        }

        const variants = await varResp.json();
        const metrics  = await metResp.json();
        const manifest = await manResp.json();
        console.log("Datos cargados:", { variants: variants.length, metrics: Object.keys(metrics).length, manifest: Object.keys(manifest).length });

        // ── KPIs calculados ───────────────────────────────────────────────
        const totalVariants    = variants.length;
        const pathogenicas     = variants.filter(v => v.classification && (v.classification.includes('Patog') || v.classification === 'Pathogenic')).length;
        const hotspots         = variants.filter(v => v.classification === 'Hotspot').length;
        const topGene          = variants[0];

        const totalStructures  = Object.values(manifest).filter(Boolean).length;
        const totalManifest    = Object.keys(manifest).length;

        const withVina         = Object.values(metrics).filter(m => m.vina_scores && Object.keys(m.vina_scores).length > 0);
        const dockedVariants   = withVina.length;

        // Mejor score global (más bajo es mejor)
        let bestScore = 0, bestDrug = '--', bestVariant = '--';
        withVina.forEach(m => {
            if (m.best_vina_score !== undefined && m.best_vina_score < bestScore) {
                bestScore   = m.best_vina_score;
                bestDrug    = m.best_vina_drug || '--';
            }
        });
        
        // Encontrar variante del mejor score
        const bestEntry = Object.entries(metrics).find(([, m]) => m.best_vina_score === bestScore);
        if (bestEntry) bestVariant = bestEntry[0].replace('_', ' ');

        // Variante con mayor RMSD activo (mayor impacto estructural)
        let maxRMSD = 0, maxRMSDVariant = '--';
        Object.entries(metrics).forEach(([id, m]) => {
            const r = parseFloat(m.rmsd_active);
            if (!isNaN(r) && Math.abs(r) > maxRMSD) {
                maxRMSD = Math.abs(r);
                maxRMSDVariant = id.replace('_', ' ');
            }
        });

        // Candidatos repositionables globales (Δ ≤ -0.5 en algún fármaco)
        const baseWT = { Alpelisib: -8.2, Abemaciclib: -8.3, Palbociclib: -8.1, Ribociclib: -8.0,
                         Copanlisib: -7.9, Capivasertib: -7.4, Everolimus: -7.6, Gedatolisib: -7.8,
                         Taselisib: -7.5, Pictilisib: -7.2, Buparlisib: -6.9 };
        let reposCount = 0;
        withVina.forEach(m => {
            Object.entries(m.vina_scores).forEach(([drug, score]) => {
                if (score !== null && baseWT[drug] && (score - baseWT[drug]) <= -0.5) reposCount++;
            });
        });

        // ── Render KPI cards ──────────────────────────────────────────────
        setText('db-total-variants',    totalVariants);
        setText('db-pathogenicas',      pathogenicas);
        setText('db-hotspots',          hotspots);
        if (topGene) {
            setText('db-top-gene',      `${topGene.hugo_symbol} ${topGene.protein_change}`);
            setText('db-top-freq',      `${topGene.frequency} casos · ${topGene.cohort_percentage?.toFixed(1)}%`);
        }
        setText('db-structures',        `${totalStructures} / ${totalManifest}`);
        setText('db-docked',            dockedVariants);
        setText('db-best-score',        `${bestScore} kcal/mol`);
        setText('db-best-drug',         bestDrug);
        setText('db-best-variant',      bestVariant);
        setText('db-max-rmsd',          `${maxRMSD.toFixed(2)} Å`);
        setText('db-max-rmsd-variant',  maxRMSDVariant);
        setText('db-repos',             reposCount);

        // ── Pipeline steps dinámicos ──────────────────────────────────────
        const m3done  = totalStructures >= 20;
        const m3pct   = totalManifest > 0 ? Math.round((totalStructures / totalManifest) * 100) : 0;
        const m4done  = dockedVariants >= 10;
        const m4pct   = totalVariants > 0 ? Math.round((dockedVariants / totalVariants) * 100) : 0;
        const m5done  = m3done && m4done;

        updateStep('step-m3', m3done ? 'done' : 'active', `${m3pct}%`);
        updateStep('step-m4', m4done ? 'done' : m3done ? 'active' : '', `${m4pct}%`);
        updateStep('step-m5', m5done ? 'done' : m4done ? 'active' : '', '');

        // ── Mini tabla Top-5 variantes ────────────────────────────────────
        const tbody = document.getElementById('db-top-table');
        if (tbody) {
            tbody.innerHTML = variants.slice(0, 5).map(v => {
                const m   = metrics[`${v.hugo_symbol}_${v.protein_change}`] || {};
                const vina = m.best_vina_score ? `${m.best_vina_score} kcal/mol` : '--';
                const drug = m.best_vina_drug || '--';
                const isPatho = v.classification && (v.classification.includes('Patog') || v.classification === 'Pathogenic');
                const cls  = isPatho
                    ? '<span style="color:#f87171;font-weight:600;">Patogénica</span>'
                    : `<span style="color:#fbbf24;">${v.classification || 'VUS'}</span>`;
                return `
                <tr style="border-bottom:1px solid var(--border);">
                    <td style="padding:10px 8px;font-weight:700;color:var(--accent);">#${v.rank}</td>
                    <td style="padding:10px 8px;"><span class="gene-badge">${v.hugo_symbol}</span></td>
                    <td style="padding:10px 8px;font-family:monospace;">${v.protein_change}</td>
                    <td style="padding:10px 8px;font-weight:600;">${v.frequency}</td>
                    <td style="padding:10px 8px;">${cls}</td>
                    <td style="padding:10px 8px;font-family:monospace;color:var(--success);">${vina}</td>
                    <td style="padding:10px 8px;font-size:0.8rem;color:var(--text-tertiary);">${drug}</td>
                </tr>`;
            }).join('');
        }
        console.log("Dashboard cargado exitosamente.");
    } catch(e) {
        console.error('Dashboard error:', e);
    }
};

function setText(id, val) {
    const el = document.getElementById(id);
    if (el) el.textContent = val;
}

function updateStep(id, status, pct) {
    const el = document.getElementById(id);
    if (!el) return;
    el.className = `step ${status}`;
    const pctEl = el.querySelector('.step-pct');
    if (pctEl && pct) pctEl.textContent = pct;
}
