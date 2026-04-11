(function () {
    const MANIFEST_URL = 'gromacs/results/gromacs_batch_manifest.json';
    const METRICS_URL = 'gromacs/results/gromacs_metrics.json';

    const fmt = (value, digits = 2) => {
        const num = Number(value);
        if (!Number.isFinite(num)) return '--';
        return num.toLocaleString('es-MX', {
            minimumFractionDigits: digits,
            maximumFractionDigits: digits
        });
    };

    async function fetchJson(url) {
        const response = await fetch(url);
        if (!response.ok) throw new Error(`No se pudo cargar ${url}`);
        return response.json();
    }

    function normalizePairs(payload) {
        const rawPairs = payload?.pairs || [];
        if (Array.isArray(rawPairs)) return rawPairs;
        return Object.entries(rawPairs).map(([tag, row]) => ({ tag, ...(row || {}) }));
    }

    function safeEl(id) {
        return document.getElementById(id);
    }

    function renderKpis(manifest, metrics) {
        const manifestPairs = Array.isArray(manifest?.pairs) ? manifest.pairs : [];
        const metricPairs = normalizePairs(metrics);
        const completed = metricPairs.filter((row) => row.status === 'completed');
        const failed = metricPairs.filter((row) => row.status === 'error' || row.status === 'failed');
        const best = completed
            .filter((row) => Number.isFinite(Number(row.mmgbsa_delta_g)))
            .sort((a, b) => Number(a.mmgbsa_delta_g) - Number(b.mmgbsa_delta_g))[0];

        if (safeEl('gmx-kpi-variants')) safeEl('gmx-kpi-variants').textContent = fmt(manifest?.summary?.variants_considered || 0, 0);
        if (safeEl('gmx-kpi-pairs')) safeEl('gmx-kpi-pairs').textContent = fmt(manifestPairs.length, 0);
        if (safeEl('gmx-kpi-completed')) safeEl('gmx-kpi-completed').textContent = fmt(completed.length, 0);
        if (safeEl('gmx-kpi-best')) safeEl('gmx-kpi-best').textContent = best ? `${fmt(best.mmgbsa_delta_g)} kcal/mol` : '--';
        if (safeEl('gmx-kpi-best-sub')) safeEl('gmx-kpi-best-sub').textContent = best ? `${best.variant_id} / ${best.drug_name}` : 'Sin MM-GBSA todavia';
        if (safeEl('gmx-status-note')) {
            safeEl('gmx-status-note').textContent = `Completados: ${completed.length} · Fallidos: ${failed.length} · Politica actual: top-${manifest?.selection_policy?.top_n_per_variant || 5} por variante.`;
        }
    }

    function renderManifestTable(manifest) {
        const body = safeEl('gmx-manifest-body');
        if (!body) return;
        const rows = Array.isArray(manifest?.pairs) ? manifest.pairs : [];
        body.innerHTML = '';

        rows.slice(0, 30).forEach((row, index) => {
            const inputState = row.available_inputs
                ? Object.entries(row.available_inputs).filter(([, ok]) => ok).length
                : 0;
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td style="padding:9px 8px;color:var(--text-tertiary);">${index + 1}</td>
                <td style="padding:9px 8px;"><code>${row.variant_id}</code></td>
                <td style="padding:9px 8px;">${row.drug_name}</td>
                <td style="padding:9px 8px;font-weight:700;color:var(--success);">${fmt(row.vina_score_kcalmol)}</td>
                <td style="padding:9px 8px;">${inputState}/7</td>
                <td style="padding:9px 8px;"><code>${row.recommended_colab_run_dir || '--'}</code></td>
            `;
            body.appendChild(tr);
        });
    }

    function renderMetricsTable(metrics) {
        const body = safeEl('gmx-metrics-body');
        if (!body) return;
        const rows = normalizePairs(metrics).sort((a, b) => {
            const av = Number.isFinite(Number(a.mmgbsa_delta_g)) ? Number(a.mmgbsa_delta_g) : 999;
            const bv = Number.isFinite(Number(b.mmgbsa_delta_g)) ? Number(b.mmgbsa_delta_g) : 999;
            return av - bv;
        });
        body.innerHTML = '';

        if (!rows.length) {
            body.innerHTML = '<tr><td colspan="8" style="padding:18px;color:var(--text-tertiary);text-align:center;">Todavia no hay resultados MD. Ejecuta la libreta de Colab y vuelve a ingerir el JSON.</td></tr>';
            return;
        }

        rows.forEach((row) => {
            const statusColor = row.status === 'completed'
                ? 'var(--success)'
                : (row.status === 'error' || row.status === 'failed' ? 'var(--danger)' : 'var(--warning)');
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td style="padding:9px 8px;"><code>${row.variant_id || '--'}</code></td>
                <td style="padding:9px 8px;">${row.drug_name || '--'}</td>
                <td style="padding:9px 8px;">${fmt(row.vina_score_kcalmol)}</td>
                <td style="padding:9px 8px;font-weight:700;color:${statusColor};">${row.mmgbsa_delta_g != null ? fmt(row.mmgbsa_delta_g) : '--'}</td>
                <td style="padding:9px 8px;">${row.rmsd_mean_nm != null ? fmt(row.rmsd_mean_nm, 3) : '--'}</td>
                <td style="padding:9px 8px;">${row.rmsf_mean_nm != null ? fmt(row.rmsf_mean_nm, 3) : '--'}</td>
                <td style="padding:9px 8px;">${row.hbonds_mean != null ? fmt(row.hbonds_mean, 2) : '--'}</td>
                <td style="padding:9px 8px;color:${statusColor};font-weight:700;">${row.status || 'unknown'}</td>
            `;
            body.appendChild(tr);
        });
    }

    function renderVariantBoard(metrics) {
        const container = safeEl('gmx-variant-board');
        if (!container) return;
        const rows = normalizePairs(metrics);
        const byVariant = {};
        rows.forEach((row) => {
            if (!row.variant_id) return;
            byVariant[row.variant_id] = byVariant[row.variant_id] || [];
            byVariant[row.variant_id].push(row);
        });

        const cards = Object.entries(byVariant).map(([variantId, items]) => {
            const best = items
                .filter((row) => Number.isFinite(Number(row.mmgbsa_delta_g)))
                .sort((a, b) => Number(a.mmgbsa_delta_g) - Number(b.mmgbsa_delta_g))[0];
            const completed = items.filter((row) => row.status === 'completed').length;
            return { variantId, best, completed, total: items.length };
        }).sort((a, b) => (a.best?.mmgbsa_delta_g ?? 999) - (b.best?.mmgbsa_delta_g ?? 999));

        container.innerHTML = '';
        if (!cards.length) {
            container.innerHTML = '<div style="color:var(--text-tertiary);font-size:0.86rem;">No hay pares con dinamica molecular ingeridos todavia.</div>';
            return;
        }

        cards.forEach((card) => {
            const el = document.createElement('div');
            el.className = 'card';
            el.style.cssText = 'padding:16px;background:rgba(255,255,255,0.02);';
            el.innerHTML = `
                <div style="display:flex;justify-content:space-between;gap:12px;align-items:center;margin-bottom:8px;">
                    <strong>${card.variantId}</strong>
                    <span style="font-size:0.75rem;color:var(--text-tertiary);">${card.completed}/${card.total} completos</span>
                </div>
                <div style="font-size:0.82rem;color:var(--text-secondary);line-height:1.6;">
                    Mejor par dinamico: <strong style="color:var(--text-primary);">${card.best ? card.best.drug_name : 'Pendiente'}</strong><br>
                    MM-GBSA: <strong style="color:var(--success);">${card.best ? `${fmt(card.best.mmgbsa_delta_g)} kcal/mol` : '--'}</strong><br>
                    RMSD ligando: ${card.best?.rmsd_mean_nm != null ? `${fmt(card.best.rmsd_mean_nm, 3)} nm` : '--'}
                </div>
            `;
            container.appendChild(el);
        });
    }

    function renderInstructions(manifest, metrics) {
        const box = safeEl('gmx-instructions');
        if (!box) return;
        const rows = normalizePairs(metrics);
        const hasResults = rows.length > 0;
        box.innerHTML = `
            <ol style="margin:0;padding-left:18px;display:grid;gap:8px;color:var(--text-secondary);font-size:0.84rem;line-height:1.6;">
                <li>Ejecuta <code>python gromacs/build_gromacs_batch.py --top_n_per_variant ${manifest?.selection_policy?.top_n_per_variant || 5}</code> para regenerar el manifiesto.</li>
                <li>Abre <code>gromacs/01_prepare_gromacs_batch_from_vina.ipynb</code> y <code>gromacs/02_run_gromacs_md_colab.ipynb</code> en Colab desde Google Drive.</li>
                <li>Guarda el JSON final en <code>gromacs/results/gromacs_metrics.json</code> dentro del repo sincronizado.</li>
                <li>Corre <code>python gromacs/ingest_gromacs_results.py</code> en local para actualizar <code>data/structural_metrics.json</code>.</li>
                <li>Recarga la app. Estado actual: <strong style="color:var(--text-primary);">${hasResults ? 'hay resultados ingeridos' : 'solo hay manifiesto y placeholder'}</strong>.</li>
            </ol>
        `;
    }

    function renderError(message) {
        const el = safeEl('gmx-error');
        if (!el) return;
        el.style.display = 'block';
        el.textContent = message;
    }

    async function initGromacsView() {
        try {
            const [manifest, metrics] = await Promise.all([
                fetchJson(MANIFEST_URL),
                fetchJson(METRICS_URL)
            ]);
            const errorBox = safeEl('gmx-error');
            if (errorBox) errorBox.style.display = 'none';
            renderKpis(manifest, metrics);
            renderManifestTable(manifest);
            renderMetricsTable(metrics);
            renderVariantBoard(metrics);
            renderInstructions(manifest, metrics);
        } catch (error) {
            console.error('Error loading GROMACS view:', error);
            renderError(`No se pudo cargar la vista GROMACS: ${error.message}`);
        }
    }

    window.initGromacsView = initGromacsView;
})();
