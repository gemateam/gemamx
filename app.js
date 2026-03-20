document.addEventListener('DOMContentLoaded', () => {
    window.selectedVariantId = null; // Initialization at the very beginning
    // Navigation logic
    const navItems = document.querySelectorAll('.nav-item');
    const views = document.querySelectorAll('.view');

    navItems.forEach(item => {
        item.addEventListener('click', () => {
            const viewId = item.getAttribute('data-view');
            if (!viewId) return;

            // Update active nav
            navItems.forEach(i => i.classList.remove('active'));
            views.forEach(v => v.classList.remove('active'));
            
            item.classList.add('active');
            const targetView = document.getElementById(`view-${viewId}`);
            if (targetView) targetView.classList.add('active');

            if (viewId === 'variants') {
                loadData();
            } else if (viewId === 'raw-data') {
                currentIndex = 0; // Reset pagination
                loadRawData();
            } else if (viewId === 'analysis') {
                if (window.loadAnalysisData) window.loadAnalysisData();
            } else if (viewId === 'structure') {
                if (window.initStructuralView) window.initStructuralView();
            } else if (viewId === 'docking') {
                if (window.initDockingView) window.initDockingView();
            } else if (viewId === 'report') {
                if (window.initReportView) window.initReportView();
            }
        });
    });

    // Pagination for Raw Data
    let rawDataStore = [];
    let currentIndex = 0;
    const PAGE_SIZE = 50;

    async function loadRawData() {
        const body = document.getElementById('raw-data-body');
        const loadMoreBtn = document.getElementById('btn-load-more');
        
        if (currentIndex === 0) body.innerHTML = '';
        
        try {
            if (rawDataStore.length === 0) {
                const response = await fetch('data/raw_mutations.json');
                rawDataStore = await response.json();
            }

            const chunk = rawDataStore.slice(currentIndex, currentIndex + PAGE_SIZE);
            chunk.forEach(item => {
                const row = document.createElement('tr');
                row.innerHTML = `
                    <td><span class="gene-badge">${item.hugo_symbol}</span></td>
                    <td>${item.variant_classification}</td>
                    <td><code>${item.protein_change}</code></td>
                    <td><span style="font-size:0.8rem; color:var(--text-tertiary)">${item.sample_id}</span></td>
                    <td>${item.chromosome}</td>
                    <td>${item.start_position}</td>
                `;
                body.appendChild(row);
            });

            currentIndex += PAGE_SIZE;
            
            if (currentIndex >= rawDataStore.length) {
                loadMoreBtn.style.display = 'none';
            } else {
                loadMoreBtn.style.display = 'inline-block';
            }
        } catch (err) {
            console.error("Error loading raw data:", err);
            body.innerHTML = `<tr><td colspan="6" style="text-align:center; color:var(--accent)">Error al cargar dataset: ${err.message}</td></tr>`;
        }
    }

    document.getElementById('btn-load-more').addEventListener('click', loadRawData);

    // Data Fetching and Filtering Logic
    const tableBody = document.getElementById('variants-body');
    const searchGene = document.getElementById('search-gene');
    const searchVariant = document.getElementById('search-variant');
    const filterClass = document.getElementById('filter-class');
    const variantCountBadge = document.getElementById('variant-count-badge');
    
    let allVariants = [];

    async function loadData() {
        try {
            const response = await fetch('data/top_variants.json');
            allVariants = await response.json();
            renderVariants(allVariants);
        } catch (error) {
            console.error('Error loading variants data:', error);
        }
    }

    function renderVariants(data) {
        if (!tableBody) return;
        tableBody.innerHTML = '';
        
        if (variantCountBadge) {
            variantCountBadge.textContent = `MOSTRANDO: ${data.length}`;
        }

        if (data.length === 0) {
            tableBody.innerHTML = '<tr><td colspan="8" style="text-align:center; padding: 40px; color: var(--text-tertiary);">No se encontraron variantes con los filtros seleccionados.</td></tr>';
            return;
        }

        const maxFreq = allVariants.length > 0 ? allVariants[0].frequency : 100;

        data.forEach(item => {
            const row = document.createElement('tr');
            row.innerHTML = `
                <td style="color: var(--text-tertiary)">${item.rank}</td>
                <td><span class="gene-badge">${item.hugo_symbol}</span></td>
                <td style="font-family: monospace; font-weight: 500;">${item.protein_change}</td>
                <td>
                    <div style="display: flex; align-items: center; gap: 10px;">
                        <div class="frequency-bar-bg">
                            <div class="frequency-bar-fill" style="width: ${(item.frequency / maxFreq) * 100}%"></div>
                        </div>
                        <span style="font-size: 0.8rem; font-weight: 600;">${item.frequency}</span>
                    </div>
                </td>
                <td style="font-weight: 600; color: #fff;">${item.cohort_percentage ? item.cohort_percentage.toFixed(1) + '%' : '0%'}</td>
                <td style="font-size: 0.85rem; color: var(--text-secondary);">Chr${item.chromosome}:${item.start_position}</td>
                <td>
                    <span style="font-size: 0.75rem; padding: 4px 8px; border-radius: 4px; background: ${item.classification === 'VUS' ? 'rgba(245, 158, 11, 0.1)' : 'rgba(239, 68, 68, 0.1)'}; color: ${item.classification === 'VUS' ? '#fbbf24' : '#f87171'}; font-weight: 600;">
                        ${item.classification || 'VUS'}
                    </span>
                </td>
                <td>
                    <div style="display: flex; gap: 5px;">
                        <a href="data/fasta/${item.hugo_symbol}_WT.fasta" target="_blank" class="fasta-btn">WT</a>
                        <a href="data/fasta/${item.hugo_symbol}_${item.protein_change}.fasta" target="_blank" class="fasta-btn" style="background: var(--accent); color: #fff;">MUT</a>
                    </div>
                </td>
            `;
            row.addEventListener('click', () => {
                window.selectedVariantId = `${item.hugo_symbol}_${item.protein_change}`;
                document.querySelectorAll('#variants-body tr').forEach(r => r.style.background = '');
                row.style.background = 'rgba(13, 115, 119, 0.1)';
            });
            
            tableBody.appendChild(row);
        });

        // Select first row of the filtered set
        if (tableBody.firstChild && data.length > 0) {
            tableBody.firstChild.click();
        }
    }

    function applyFilters() {
        const geneQuery = searchGene ? searchGene.value.toLowerCase().trim() : '';
        const variantQuery = searchVariant ? searchVariant.value.toLowerCase().trim() : '';
        const classVal = filterClass ? filterClass.value : 'ALL';

        const filtered = allVariants.filter(v => {
            const matchesGene = !geneQuery || v.hugo_symbol.toLowerCase().includes(geneQuery);
            const matchesVariant = !variantQuery || v.protein_change.toLowerCase().includes(variantQuery);
            const matchesClass = (classVal === 'ALL') || (v.classification === classVal);
            return matchesGene && matchesVariant && matchesClass;
        });

        renderVariants(filtered);
    }

    if (searchGene) searchGene.addEventListener('input', applyFilters);
    if (searchVariant) searchVariant.addEventListener('input', applyFilters);
    if (filterClass) filterClass.addEventListener('change', applyFilters);

    loadData();

    // Sync button interaction
    const btnSync = document.getElementById('btn-sync');
    if (btnSync) {
        btnSync.addEventListener('click', () => {
            btnSync.innerHTML = '<span>Sincronizando...</span>';
            setTimeout(() => {
                btnSync.innerHTML = '<span>Sincronizado ✓</span>';
                btnSync.style.background = 'var(--success)';
            }, 1500);
        });
    }

    // Analyze Structure Button
    const btnAnalyze = document.getElementById('btn-analyze-structure');
    if (btnAnalyze) {
        btnAnalyze.addEventListener('click', () => {
            const variantToLoad = window.selectedVariantId;
            showView('structure');
            
            // Transfer selection to structural view dropdown
            if (variantToLoad) {
                setTimeout(() => {
                    const select = document.getElementById('struct-variant-select');
                    if (select) {
                        select.value = variantToLoad;
                        if (select.value === variantToLoad) {
                            select.dispatchEvent(new Event('change'));
                        }
                    }
                }, 300);
            }
        });
    }

    function showView(viewId) {
        const navItem = document.querySelector(`.nav-item[data-view="${viewId}"]`);
        if (navItem) navItem.click();
    }

    // Continue to Docking Button
    const btnNextDocking = document.getElementById('btn-next-docking');
    if (btnNextDocking) {
        btnNextDocking.onclick = () => showView('docking');
    }


    // --- 3D Viewer Navigation Helpers ---
    window.zoomViewer = (viewerKey, ratio) => {
        const plugin = { wt: window.wtViewer, mut: window.mutViewer, docking: window.dockingViewer }[viewerKey];
        if (!plugin) return;

        // PDBe Molstar v3 plugin instance (the base Mol* plugin)
        const molstar = plugin.plugin;
        if (!molstar?.canvas3d) return;

        const camera = molstar.canvas3d.camera;
        const current = camera.state;
        
        // ratio > 1 = zoom in (reduce radius), ratio < 1 = zoom out (increase radius)
        const newRadius = current.radius * (ratio > 1 ? 0.75 : 1.35);

        molstar.canvas3d.requestCameraReset({
            snapshot: {
                ...current,
                radius: Math.max(newRadius, 1)
            },
            durationMs: 250
        });
    };

    window.resetViewer = (viewerKey) => {
        const plugin = { wt: window.wtViewer, mut: window.mutViewer, docking: window.dockingViewer }[viewerKey];
        const molstar = plugin?.plugin;
        if (molstar?.canvas3d) {
            molstar.canvas3d.requestCameraReset({ durationMs: 400 });
        }
    };

    window.toggleExpand = (viewerKey) => {
        const panelId = viewerKey === 'docking' ? 'docking-viewer-container' : `${viewerKey}-viewer-panel`;
        const panel = document.getElementById(panelId);
        if (!panel) return;

        const isExpanded = panel.classList.contains('force-expanded');
        panel.classList.toggle('force-expanded', !isExpanded);
        
        // Ensure WebGL canvas re-fills correctly
        requestAnimationFrame(() => {
            requestAnimationFrame(() => {
                window.dispatchEvent(new Event('resize'));
            });
        });
    };
});
