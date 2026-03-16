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
            item.classList.add('active');

            // Update visible view
            views.forEach(v => {
                v.classList.remove('active');
                if (v.id === `view-${viewId}`) {
                    v.classList.add('active');
                }
            });
        });
    });

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
                        <td>
                            <span style="font-size: 0.75rem; padding: 4px 8px; border-radius: 4px; background: ${relevance === 'VUS' ? 'rgba(245, 158, 11, 0.1)' : 'rgba(239, 68, 68, 0.1)'}; color: ${relevance === 'VUS' ? '#fbbf24' : '#f87171'}; font-weight: 600;">
                                ${relevance}
                            </span>
                        </td>
                    `;
                    tableBody.appendChild(row);
                });
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
});
