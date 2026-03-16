document.addEventListener('DOMContentLoaded', () => {
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

    // Data Fetching Logic
    const tableBody = document.getElementById('variants-body');
    
    async function loadData() {
        try {
            const response = await fetch('data/top_variants.json');
            const variantsData = await response.json();
            
            if (tableBody) {
                tableBody.innerHTML = '';
                variantsData.forEach(item => {
                    const row = document.createElement('tr');
                    // Frequency normalization based on max (128 for PIK3CA)
                    const maxFreq = variantsData[0].frequency;
                    const relevance = item.frequency > 50 ? 'High Pathogenic' : (item.frequency > 10 ? 'Pathogenic' : 'VUS');
                    
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
                        // Visual feedback
                        document.querySelectorAll('#variants-body tr').forEach(r => r.style.background = '');
                        row.style.background = 'rgba(13, 115, 119, 0.1)';
                    });
                    
                    tableBody.appendChild(row);
                });
                
                // Select first row by default
                if (tableBody.firstChild) {
                    tableBody.firstChild.click();
                }
            }
        } catch (error) {
            console.error('Error loading variants data:', error);
        }
    }

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
            showView('structure');
        });
    }

    function showView(viewId) {
        const navItem = document.querySelector(`.nav-item[data-view="${viewId}"]`);
        if (navItem) navItem.click();
    }

    window.selectedVariantId = null;
});
